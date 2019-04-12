
pitchofav = 5*pi/180;
betaofav = 5 *pi/180;
TAS = 5;


Vf2 = [cos(pitchofav); 
       sin(pitchofav)]*TAS;

Vf = [Vf2(1)*cos(betaofav)
      Vf2(1)*sin(betaofav)
      Vf2(2)];
  
omegasRPM= (2.5:1:6.5)*1e3;
% omegasRPM= 12*1e3;
Nb = 2;
Rcm = 14;
R = 14/100; 



theta = @(r)((0.0009159*(r*Rcm).^5 - 0.04202*(r*Rcm).^4 + 0.742*(r*Rcm).^3 -...
             6.151*(r*Rcm).^2 + 21.24*(r*Rcm) + 0.1216 )*(pi/180));

C = @(r)((-0.001067*(r*Rcm).^3 - 0.02944*(r*Rcm).^2 + 0.6259.*(r*Rcm) + 0.5392)/100);



CL = @(alpha)interp1([-180 -12   -10    -5  10 12 21 180]*pi/180,...
                     [   0   0 -0.25 -0.25   1  1  0   0],alpha);
                 
CD = @(alpha)interp1([-180  -70    -10    10   30   70  180]*pi/180,...
                     [ 1.3  1.3   0.05  0.05  0.7  1.3  1.3],alpha);



                 
y = linspace(2,14,100);
alpha1 = linspace(-40,40,100);
alpha2 = linspace(-70,70,100);

figure(1);
clf
ax = [subplot(2,2,1);
      subplot(2,2,2);
      subplot(2,2,3);
      subplot(2,2,4)];
  
plot(ax(1),y,theta(y/max(y))*180/pi);
title(ax(1),'pitch');
xlim(ax(1),[min(y) max(y)])
ylim(ax(1),[5 30])

plot(ax(2),y,C(y/max(y))*100);
title(ax(2),'chord');
xlim(ax(2),[min(y) max(y)])
ylim(ax(2),[0.4 3.75])

plot(ax(3),alpha1,CL(alpha1*pi/180));
title(ax(3),'CL');
xlim(ax(3),[min(alpha1) max(alpha1)])

plot(ax(4),alpha2,CD(alpha2*pi/180));
title(ax(4),'CD');
xlim(ax(4),[min(alpha2) max(alpha2)])





% [T, a, P, rho] = atmosisa(0);

rho = 1.225;
Nw =  length(omegasRPM);


FMs = zeros(6,Nw);
vis = zeros(Nw,1);

% phis = zeros(Nw,1);
% vrs = zeros(Nw,1);
% rs = zeros(Nw,1);
% psis = zeros(Nw,1);

figure(2);
clf;
ax2 = arrayfun(@(x)subplot(Nw,1,x),1:Nw);

figure(3);
clf;
ax3 = arrayfun(@(x)subplot(Nw,1,x),1:Nw);

for i = 1:Nw
    omega = omegasRPM(i)*2*pi/60; % rad/s
    [CFM,vis(i),vidrdp,phis,vrs,rs,psis] = HBEMT(Vf,omega,R,Nb,CL,CD,theta,C);
    
%     fac = (0.5*0.25);
    fac = 1;
    q = fac*rho*pi*R^4*omega^2;
    FMs(:,i) = CFM*q;
    
    if numel(phis) == max(size(phis))
        alphas = theta(rs)-phis';
    else
        alphas = theta(rs)-phis';
    end
    
    
    plot(ax2(i),rs,alphas*180/pi,'r',rs,theta(rs)*180/pi,'b',rs,phis*180/pi,'k');
    plot(ax3(i),rs,vidrdp,'k');
    
end

%%
fh = figure(4);
clf;
ax = subplot(1,1,1);
plot(ax(1),omegasRPM/1e3,FMs(1,:));
ax.YTick = 0:2:8;
ax.YLim = [0 8];
ax.XTick = 2:1:7;
ax.XLim = [2 7];

grid(ax,'on');









