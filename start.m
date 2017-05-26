% INITIALIZATION AND START UP
clc; close all; clear;

global g deg air loc bm act wind X 

mis = [  50,   0,   0,  10;   
        100,   0,   0,   0;
        100,  50,  pi,  10;
          0,  60,  pi,  10;
          0,  65,  pi,  50;]; % (Y, X, heading angle, T[s]) format

%initial conditions for DRONI 1
X7_ini1=[40 0 -77 1 0 0 0]'; %initial values for [N E D e0 e1 e2 e3];
X6_ini1=[2 0   0 0 0 0]'; %initial values for speed [u v w p q r];

%initial conditions for DRONI 2
X7_ini2=[30 -10 -77 1 0 0 0]'; %initial values for [N E D e0 e1 e2 e3];
X6_ini2=[  0   0   0 0 0 0]'; %initial values for speed [u v w p q r];

%initial conditions for DRONI 3
X7_ini3=[20 10 -77 1 0 0 0]'; %initial values for [N E D e0 e1 e2 e3];
X6_ini3=[  0  0   0 0 0 0]'; %initial values for speed [u v w p q r];

angFormation1 = 0;
angFormation2 = deg2rad(-45);
angFormation3 = deg2rad(45);

[loc,air,bm,act,wind,X]=db_droni('bm65_DRONI.xls'); % airship parameters

wind.ws=3; % wind speed
wind.wh=30*deg; % wind incidence angle (from North)
wind.turb=1; %flag to activate turbulence
wind.sigg=1; %turbulence std deviation = 1 (max. sigg=7 -> storm)

Di=-loc.hG-50; %50 m height

% mission parameters for waypoint trajectory
ui=9;
Di=-77;
uref=ui;

hsim=-Di;   du1=.5;  dh1=1;
u2=5;    h2=5;      du2=.125;   dh2=.5;
u3=5;    h3=2;      du3=.1;     dh3=.2;

usim=ui;du1=.05*usim;dh1=.1*usim;
Ei=0;
Ni=0;

%mission definition
%
% (time is already introduced but not used yet, should be necessary for hovering)

% mis=[N E u h du dh tt t]             %%%%% u is Vt
%
% .N,E,h are way points coordinates (beginning of segment)(in meters from initial point)
% .u is speed for segment beginning at way point :is air speed by now (Vt) in m/s
% .du,dh are lon. reference rates for transitions
% .tt is angular rotation for circular segment (0 is straight, pi is half circle)
% .t is time (not used by now)
%---------------------------------------
% mis = [0  0 usim hsim du1 dh1  0 0;
%      100  0 usim hsim du1 dh1 pi 0;
%      100 50 usim hsim du1 dh1  0 0;
%        0 50 usim hsim du1 dh1 pi 0;
%        0  0 usim hsim du1 dh1  0 0;];

taugui = 4.5;
tau    = taugui;
X.Ts   = 1/50;
Ts     = 0.1;
dt0    = 0.1;
xmin = -150; xmax = 115; ymin = -120; ymax = 105; st = -1;

% -------------------------
% Calculating controller for DRONI 2 & 3
% -------------------------
kp = 1; ka = 1; kb = 1;
A = [-kp 0 0; 0 -(ka-kp) -kb; 0 -kp 0];
B = [-1 0; 0 -1; 0 0];
Q = [1 0 0; 0 1 0; 0 0 1];
R = [1 0; 0 1];

Kg = lqr(A,B,Q,R);
closedSys = A-B*Kg;
a  = 2;
% kp = -closedSys(1,1)/a;       % LQR
% kb = -closedSys(2,3)/a*1.5;   % LQR
% ka = kp - closedSys(2,2)/a;   % LQR
kp = 0.5;
ka = 5;
kb = (2/pi*kp - ka + 0.1)*3/5;

% -------------------------
% Stability check
% -------------------------
if kp > 0 && kb < 0 && ka-kp > 0
    stability = 1;
else
    stability = 0;
end

% -------------------------
% Robust stability check
% -------------------------
if kp > 0 && kb < 0 && ka + 5/3*kb - 2/pi*kp > 0
    robustStability = 1;
else
    robustStability = 0;
end

sim('sbm65_kin'); % call simulink model

% Saving data
out.time         = speedData.time;
out.controllerU  = speedData.signals(1).values;
out.controllerV  = speedData.signals(2).values;
out.controllerR  = speedData.signals(3).values;
out.controllerVt = speedData.signals(4).values;
out.realU        = linearSpeedData.signals(1).values;
out.realV        = linearSpeedData.signals(2).values;
out.realR        = angularSpeedData.signals(3).values;
out.realVt       = speedData.signals(4).values;
out.trajectories = trajectories;
% out.distance     = distanceTargetData.signals.values(:);
out.beta         = airDynamicsData.signals.values;

plotInfo = plotAirships(2, 0:30:150, trajectories(:,1:9), wind, out.time);
hold on
% plotInfo(4) = plot(out.trajectories(:,11),out.trajectories(:,10),'black-.'); % Plot moving target
plotInfo(4) = plot(mis(:,2),mis(:,1),'x', 'color', 'black'); % Plot mis points
for misLine = 1:size(mis,1) % Plot hover points 
    if mis(misLine,4) ~= 0
        viscircles([mis(misLine,2),mis(misLine,1)],5) % ([x,y],radius)
    end
end
plotInfo(5) = plot([X7_ini1(2);X7_ini2(2);X7_ini3(2)],[X7_ini1(1);X7_ini2(1);X7_ini3(1)],'x', 'color', 'black');
hold off

axisSet = [get(gca,'xlim'),get(gca,'ylim')]; % For the movie
legend(plotInfo, 'Droni 1', 'Droni 2', 'Droni 3', 'Ground target','Start Points');

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3.76*1.4 2.82*1.4])
print('waypointHover','-dpng','-r100')

figure(3)
subplot(2,1,1)
plot(out.time,out.controllerU,out.time,out.realU)
ylabel('$$u_a$$ [m/s]','Interpreter','latex','FontSize',13)
xlabel('Simulation time [s]')
legend('Controller output','Actual velocity')

subplot(2,1,2)
plot(out.time,out.controllerR,out.time,out.realR)
ylabel('$$r$$ [m/s]','Interpreter','latex','FontSize',13)
xlabel('Simulation time [s]')

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3.76*1.4 2.82*1.4])
print('waypointHoverSpeed','-dpng','-r100')

figure(4)
% subplot(2,1,1)
subplot(3,1,1)
plot(out.time,out.distance, out.time,ones(size(out.time))*5)
ylabel({'Distance to',' target [m]'},'Interpreter','latex','FontSize',13)
% xlabel('Simulation time [s]')

% subplot(2,2,3)
subplot(3,1,2)
plot(out.time,out.realVt)
ylabel('$$V_t$$ [m/s]','Interpreter','latex','FontSize',13)
% xlabel('Simulation time [s]')

% subplot(2,2,4)
subplot(3,1,3)
plot(out.time,out.beta)
ylabel('$$\beta$$ [m/s]','Interpreter','latex','FontSize',13)
% xlabel('Simulation time [s]')

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3.76*1.4 2.82*1.4])
print('movingTargetDynamics','-dpng','-r100')

% axisSet = [get(gca,'xlim'),get(gca,'ylim')];
out.movieFrames = plotMovie(5,1:0.5:150,out.trajectories,wind,out.time,axisSet);