function varargout = HBEMT(Vf,omega,R,Nb,CL,CD,theta,C,varargin)
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

p = inputParser;

% radial position [m]
p.addParameter('r0'   ,2/100       ,@isscalar);
p.addParameter('rf'   ,14/100      ,@isscalar);
% blade azimuth in rotor-wind axis (V axis - H axis plane along Vxy-fig.5) [rad]
p.addParameter('psi0' ,0           ,@isscalar);
p.addParameter('psif' ,2*pi        ,@isscalar);
% Integration Mode [-]
p.addParameter('integrationmode',3 ,@isscalar);
% Number of places to evaluate integral (radial and azimuth axis) [-]
p.addParameter('Nr'  ,-1           ,@isscalar);
p.addParameter('Naz' ,-1           ,@isscalar);

p.parse(varargin{:});

r0              = p.Results.r0/R; % radial position as ratio of total radius (R)[-]
rf              = p.Results.rf/R; % radial position as ratio of total radius (R)[-]
psi0            = p.Results.psi0;
psif            = p.Results.psif;
integrationmode = p.Results.integrationmode;
Nr              = p.Results.Nr;
Naz             = p.Results.Naz;

switch integrationmode
    case {1,2}
    % for simpson/trapozoidal integration
        if Naz < 0
            Naz = 360;
        end
        if Nr < 0
            Nr = 10;
        end
        r = linspace(r0,rf,Nr);             % radial position as ratio of total radius (R)[-]
        
        if (2*pi - (psif - psi0)) <= eps('single')
            psi0 = psi0 + abs(psif - psi0)/Naz;
        end
        
        psi = linspace(psi0,psif,Naz); % blade azimuth in rotor-wind axis i.e. V axis along Vxy  - fig 5
    case 3
    % gauss legendre integration
        if Naz < 0
            Naz = 10;
        end
        if Nr < 0
            Nr = 6;
        end
    
        [r,rc]=lgwt(Nr,r0,rf);
        r= r';
        Nr  = length(r);
        
        psif = 2*pi;
        psi0  = 0;
        [psi,psic]=lgwt(Naz,psi0,psif);
        psi = psi';
        Naz = length(psi);
end

Vxy = sqrt(sum(Vf(1:2).^2));
Vz = Vf(3);

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
                            

        lambdai(i,j) = brentDekker(aux,0.01,1.0);
        [CTdrdp(i,j),CVdrdp(i,j),CHdrdp(i,j),...
         CMTdrdp(i,j),CMVdrdp(i,j),CMHdrdp(i,j),...
         phidrdp(i,j),vrdrdp(i,j)] ...
            = bladeElementForcesAndMoments(lambdai(i,j),r(i),psi(j),...
                lambdac,mu,sigmaPrime,CL,CD,theta,C);
    end
end


vidrdp = lambdai*omega*R;

switch integrationmode
    case 1
        % trapazoidal intergatrion
        vi =   trapz(psi,trapz(r,vidrdp));
        CFM = [trapz(psi,trapz(r,CTdrdp));
               trapz(psi,trapz(r,CVdrdp));
               trapz(psi,trapz(r,CHdrdp));
               trapz(psi,trapz(r,CMTdrdp));
               trapz(psi,trapz(r,CMVdrdp));
               trapz(psi,trapz(r,CMHdrdp))];
    case 2
        % simpson intergatrion
        dpsi = diff(psi(1:2));
        dr = diff(r(1:2));
        vi =   simpson(dpsi,simpson(dr,vidrdp));
        CFM = [simpson(dpsi,simpson(dr,CTdrdp));
               simpson(dpsi,simpson(dr,CVdrdp));
               simpson(dpsi,simpson(dr,CHdrdp));
               simpson(dpsi,simpson(dr,CMTdrdp));
               simpson(dpsi,simpson(dr,CMVdrdp));
               simpson(dpsi,simpson(dr,CMHdrdp))];
    case 3
        vi =   gaussLegendre(psic,gaussLegendre(rc,vidrdp));
        CFM = [gaussLegendre(psic,gaussLegendre(rc,CTdrdp));
               gaussLegendre(psic,gaussLegendre(rc,CVdrdp));
               gaussLegendre(psic,gaussLegendre(rc,CHdrdp));
               gaussLegendre(psic,gaussLegendre(rc,CMTdrdp));
               gaussLegendre(psic,gaussLegendre(rc,CMVdrdp));
               gaussLegendre(psic,gaussLegendre(rc,CMHdrdp))];
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


% rotor disk axis 
% X axis along rotor shaft
% Y axis is in direction of wind in plane of rotor disk 
% Z axis is lateral to wind in plane of rotor disk 
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


function I = simpson(h,y)
[n,m] = size(y);
if n>1 && m>1
    % for matrix inegrate along column and return row vector of intgrated
    % values
    I = h/3*(y(1,:)+2*sum(y(3:2:end-2,:))+4*sum(y(2:2:end,:))+y(end,:)); 
else
    I = h/3*(y(1)+2*sum(y(3:2:end-2))+4*sum(y(2:2:end))+y(end));
end
end


function I = gaussLegendre(c,y)
[n,m] = size(y);
I = zeros(1,m);

c = c(:);% ensure c is a column vector

gl = @(f)(sum(c.*f));
if n>1 && m>1
    % for matrix inegrate along column and return row vector of intgrated
    % values
    for i = 1:m
        I(i) = gl(y(:,i));
    end
else
    y = y(:);% ensure y is a column vector
    I = gl(y);
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

