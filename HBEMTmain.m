climbing = true;
pitchofav = 45*pi/180;
betaofav = 0 *pi/180;
TAS = 50;

if climbing
    Vf = [0; 0; TAS];
else
    % otherwise assume holding altitude and traveling forward
    Vf2 = [cos(pitchofav);
           sin(pitchofav)]*TAS;
       
    Vf = [Vf2(1)*cos(betaofav)
          Vf2(1)*sin(betaofav)
          Vf2(2)];
end


  
% omegasRPM= (2.5:0.5:12.5)*1e3;
omegasRPM= linspace(0,30,100)*1e3;
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





[T, a, P, rho] = atmosisa(0);

Nw =  length(omegasRPM);


FMs = zeros(6,Nw);
vis = zeros(Nw,1);

% phis = zeros(Nw,1);
% vrs = zeros(Nw,1);
% rs = zeros(Nw,1);
% psis = zeros(Nw,1);
if Nw <=4
    figure(2);
    clf;
    ax2 = arrayfun(@(x)subplot(Nw,1,x),1:Nw);
    
    figure(3);
    clf;
    ax3 = arrayfun(@(x)subplot(Nw,1,x),1:Nw);
end
tms = zeros(Nw,1);
for i = 1:Nw
    omega = omegasRPM(i)*2*pi/60; % rad/s
    tic
    [CFM,vis(i),vidrdp,phis,vrs,rs,psis] = HBEMT(Vf,omega,R,Nb,CL,CD,theta,C,...
                                                 'Nr',6,...
                                                 'Naz',10,...
                                                 'integrationmode',3);
    tms(i) = toc;
%     fac = (0.5*0.25);
    fac = 1;
    FMs(1:3,i) = CFM(1:3,1)*rho*pi*R^4*omega^2;
    FMs(4:6,i) = CFM(4:6,1)*rho*pi*R^5*omega^2;
    alphas = theta(rs)'-phis;
%     if numel(phis) == max(size(phis))
%         alphas = theta(rs)-phis';
%     else
%         alphas = theta(rs)-phis';
%     end
    
    if Nw <=4
        plot(ax2(i),rs,alphas*180/pi,'r',rs,theta(rs)*180/pi,'b',rs,phis*180/pi,'k');
        plot(ax3(i),rs,vidrdp,'k');
    end
    
end

fprintf(1,'On average HBEMT took %.8f s to run \n',mean(tms));

%%
fh = figure(5);
clf;
ax = arrayfun(@(x)subplot(3,2,x),[ 1 3 5 2 4 6]);
titleStrs = {'FT' 'FV' 'FH' 'MT' 'MV' 'MH'};



for i = 1:6
    plot(ax(i),omegasRPM/1e3,FMs(i,:));
    ax(i).Title.String = titleStrs{i};
    if i == 1 && false
        ax(i).YTick = 0:2:8;
        ax(i).YLim = [0 8];
    end
    
    ax(i).XLim = [omegasRPM(1) omegasRPM(end)]/1e3;
    ax(i).XLim = round(ax(i).XLim)+[-1 1];
    ax(i).XTick = ax(i).XLim(1):1:ax(i).XLim(end);
    grid(ax(i),'on');
end


%%
fh = figure(6);
clf
ax = axes(fh);
plot(ax,omegasRPM/1e3,vis);

d = 1;
dx = diff(omegasRPM([end-d,end])/1e3);
dy = diff(vis([end-d,end]));
a = dy/dx;
b =  vis(end) - a*omegasRPM(end)/1e3;
hold(ax,'on');
% plot(ax,omegasRPM/1e3, (a*omegasRPM/1e3) + b);
ax.Title.String = 'induced velocity [m/s]';









