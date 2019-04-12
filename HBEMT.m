function varargout = HBEMT(Vf,omega,R,Nb,CL,CD,theta,C)
% Hybrid Blade Element Momentum Theory for computing thrust in
% forward flight.
% reference: 
% Generalized Flight Dynamic Model of Quadrotor Using Hybrid Blade Element 
% Momentum Theory. DOI: 10.2514/1.C034899
% variable        
%       Vf      freestream air velocity vector
%    omega      angular speed of rotor              [rad/s]
%        R      total radius of rotor disk          [m]
%       Nb      Number of blades                    [blades]
%       CL      function w.r.t to alpha             [-]
%       CD      function w.r.t to alpha             [-]
%    theta      function w.r.t to span ratio r      [rad]
%        C      function w.r.t to span ratio r      [rad]
% freestream air velcoty vector (Vf)


Vxy = sqrt(sum(Vf(1:2).^2));
Vz = Vf(3);
if Vxy == 0 
    Naz = 1;
else
    Naz = 5;
end
Nr = 7;


r = linspace(2/14,1,Nr);            % radial position as ratio of total radius (R)[-]
psi = linspace(0,2*pi,Naz);      % blade azimuth in rotor-wind axis i.e. V axis along Vxy 




lambdac = Vz/(omega*R);
mu = Vxy/(omega*R);

sigmaPrime = Nb/(pi*R);    % multiply by chord for solidity variation



szTmp = [Nr,Naz];
lambdai = zeros(szTmp);
CTdrdp  = zeros(szTmp);
CVdrdp  = zeros(szTmp);
CHdrdp  = zeros(szTmp);
CMTdrdp = zeros(szTmp);
CMVdrdp = zeros(szTmp);
CMHdrdp = zeros(szTmp);
phidrdp = zeros(szTmp);
vrdrdp = zeros(szTmp);

% tms = [];

for i = 1:Nr
    for j = 1:Naz
        aux =@(lambdai)HBEMTaux(lambdai,r(i),psi(j),lambdac,mu,sigmaPrime,CL,...
                                CD,theta,C);
                            
%         g = 0.1;
%         while true
%             lambdai(i,j) = fzero(aux,g);
%             if ~isnan(lambdai(i,j)) && ~isinf(lambdai(i,j))
%                 break;
%             end
%             
%             fprintf(1,['At psi %f and r %f\n',...
%                        'try again. lambdai %f guess %f\n\n\n'],...
%                     psi(j)*180/pi,r(i),lambdai(i,j),g);
%             
%             if lambdai(i,j) < 0
%                 g = g*1.1;
%             else
%                 g = g*0.5;
%             end
% 
%             if g < 0.0063
%                 error('could not find lambdai');
%             end
%             
% 
%             
%         end
        lambdai(i,j) = brentDekker(aux,0.01,.2);
        [CTdrdp(i,j),CHdrdp(i,j),CVdrdp(i,j),...
         CMTdrdp(i,j),CMHdrdp(i,j),CMVdrdp(i,j),...
         phidrdp(i,j),vrdrdp(i,j)] ...
            = bladeElementForcesAndMoments(lambdai(i,j),r(i),psi(j),...
                lambdac,mu,sigmaPrime,CL,CD,theta,C);
    end
end


vidrdp = lambdai*omega*R;
if isscalar(psi)
    vi = psi*trapz(r,vidrdp);
    
    CFM = [psi*trapz(r,CTdrdp);
           psi*trapz(r,CVdrdp);
           psi*trapz(r,CHdrdp);
           psi*trapz(r,CMTdrdp);
           psi*trapz(r,CMVdrdp);
           psi*trapz(r,CMHdrdp)];
else
    vi = trapz(psi,trapz(r,vidrdp));
    
    CFM = [trapz(psi,trapz(r,CTdrdp));
           trapz(psi,trapz(r,CVdrdp));
           trapz(psi,trapz(r,CHdrdp));
           trapz(psi,trapz(r,CMTdrdp));
           trapz(psi,trapz(r,CMVdrdp));
           trapz(psi,trapz(r,CMHdrdp))];
    
end

varargout = {CFM,vi,vidrdp};

if nargout > 3
    varargout{4} = phidrdp; 
end

if nargout > 4
    varargout{5} = vrdrdp;
end

if nargout > 5
    varargout{6} = r;
end

if nargout > 6
    varargout{7} = psi;
end

end

function res = HBEMTaux(lambdai,r,psi,lambdac,mu,sigmaPrime,CL,CD,theta,C)

% calcuation momentum theory
v = lambdac+lambdai;

% calculation for blade element theory
sigma = sigmaPrime*C(r);
lambdaxyc = mu*cos(psi);
phi = atan2(v,r+lambdaxyc);
alpha = theta(r)-phi;
vr2= (r+lambdaxyc)^2+v^2;

      
res = 4*lambdai*v*r      ...                                         %   momentum theory
      -       0.5*sigma*(CL(alpha)*cos(phi)-CD(alpha)*sin(phi))*vr2; % minus   Blade Element Theory

end

function varargout = bladeElementForcesAndMoments(lambdai,r,psi,lambdac,mu,...
                        sigmaPrime,CL,CD,theta,C)

v = lambdac+lambdai;
% calculation for blade element theory
sigma = sigmaPrime*C(r);
lambdaxyc = mu*cos(psi);
phi = atan2(v,r+lambdaxyc);
alpha = theta(r)-phi;
vr2= (r+lambdaxyc)^2+v^2;

tmp = (1/(2*pi))*(1/2)*sigma*vr2;
      
varargout{1} = tmp*(CL(alpha)*cos(phi)-CD(alpha)*sin(phi));            % [CT  dr dpsi]    elemental coefficient of thrust
varargout{2} = tmp*(CL(alpha)*sin(phi)+CD(alpha)*cos(phi))*cos(psi);   % [CV  dr dpsi]    elemental coefficient of with free stream force
varargout{3} = tmp*(CL(alpha)*sin(phi)+CD(alpha)*cos(phi))*sin(psi);   % [CH  dr dpsi]    elemental coefficient of lateral to free stream force
varargout{4} = tmp*(CL(alpha)*sin(phi)+CD(alpha)*cos(phi))*r;          % [CMT dr dpsi]    elemental coefficient of yaw moment
varargout{5} = tmp*(CL(alpha)*cos(phi)-CD(alpha)*sin(phi))*r*cos(psi); % [CMV dr dpsi]    elemental coefficient of with free stream moment
varargout{6} = tmp*(CL(alpha)*cos(phi)-CD(alpha)*sin(phi))*r*sin(psi); % [CMH dr dpsi]    elemental coefficient of lateral to free stream moment

if nargout > 6
    varargout{7} = phi;
end

if nargout > 7
    varargout{8} = sqrt(vr2);
end


end

function x = brentDekker(f,a,b,varargin)
% http://people.sc.fsu.edu/~jburkardt/m_src/brent/zero.m

t = eps('double');
machep = eps('double');

sa = a;
sb = b;
fa = f ( sa );
fb = f ( sb );

c = sa;
fc = fa;
e = sb - sa;
d = e;

maxiter = 100;



for i=1:maxiter
    
    if ( abs ( fc ) < abs ( fb ) )
        
        sa = sb;
        sb = c;
        c = sa;
        fa = fb;
        fb = fc;
        fc = fa;
        
    end
    
    tol = 2.0 * machep * abs ( sb ) + t;
    m = 0.5 * ( c - sb );
    
    if ( abs ( m ) <= tol || fb == 0.0 )
        break
    end
    
    if ( abs ( e ) < tol || abs ( fa ) <= abs ( fb ) )
        
        e = m;
        d = e;
        
    else
        
        s = fb / fa;
        
        if ( sa == c )
            
            p = 2.0 * m * s;
            q = 1.0 - s;
            
        else
            
            q = fa / fc;
            r = fb / fc;
            p = s * ( 2.0 * m * q * ( q - r ) - ( sb - sa ) * ( r - 1.0 ) );
            q = ( q - 1.0 ) * ( r - 1.0 ) * ( s - 1.0 );
            
        end
        
        if ( 0.0 < p )
            q = - q;
        else
            p = - p;
        end
        
        s = e;
        e = d;
        
        if ( 2.0 * p < 3.0 * m * q - abs ( tol * q ) && p < abs ( 0.5 * s * q ) )
            d = p / q;
        else
            e = m;
            d = e;
        end
        
    end
    
    sa = sb;
    fa = fb;
    
    if ( tol < abs ( d ) )
        sb = sb + d;
    elseif ( 0.0 < m )
        sb = sb + tol;
    else
        sb = sb - tol;
    end
    
    fb = f ( sb );
    
    if ( ( 0.0 < fb && 0.0 < fc ) || ( fb <= 0.0 && fc <= 0.0 ) )
        c = sa;
        fc = fa;
        e = sb - sa;
        d = e;
    end
    
end

x = sb;

return
end

% function x = rootFind(fun,x0)
% 
% x = x0;
% fx = fun(x);
% tol = ones(size(fx))*1e-10;
% 
% while any(abs(fx) > tol)
%     dx = eps('single')*2;
%     Del_f = finDiff(fun,x,dx);
%     x = x-fx/Del_f;
%     fx = fun(x);
% end
% 
% 
% end
% 
% function Del_f = finDiff(f,X,DX)
% % Fintie Difference
% % This function computes the numerical derivative of a function f(x) using 
% % the Finite Difference method
% % If X is a vector then this computes the numerical gradiant of f(x).
% %       Del_f = [df/dx1 df/dx2 ... df/dxn]
% % if f(x) is a system of equations this computes the numerical Jacobian of
% % the system
% %      J = [df1/dx1 df1/dx2 ... df1/dxn] 
% %           df2/dx1 df2/dx2 ... df2/dxn
% %           ...                     ...
% %           dfm/dx1 dfm/dx2 ... dfm/dxn
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% if ~isvector(X)
%    error('X should be vector');
% end
% 
% N = numel(X);
% 
% if nargin < 3
%     DX = 2*sqrt(eps)*ones(size(X));
% elseif numel(DX) == 1
%     DX = DX*ones(size(X));
% elseif ~isvector(DX) && numel(DX) ~= N
%     error('the optional argument delX must be a scalar or vector of the same dimension of X');
% end
% 
% 
% f0 = f(X);
% if ~isvector(f0) || ~isscalar(f0)
%    error('f(x) should be a scalar or vector');
% end
% 
% M = numel(f0);
% Del_f = zeros(M,N);
% 
% Xp = X;
% for i = 1:N
%     Xp(i) = Xp(i)+DX(i);
%     Del_f(i) = ( f(Xp) - f0 )/DX(i);
%     Xp(i) = Xp(i)-DX(i);
% end
% 
% end