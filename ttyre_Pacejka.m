function [sx,aad,Mx,My]=ttyre_Pacejka(inp);

% function [sx,aad,Mx,My]=ttyre2(inp);
%
% define tyre 2D force parameters (and plot Fx Fy curves w/o inp)
%
% data are used in fdslip dynamic function
%091117

global Bx Cx Ex By Cy Ey Br Cr Er sxm aadm

% tyres model as in ttyre2 < Milliken / Brian Beckman
Bx=.0822;Cx=1.65;Ex=-10; %added for combined slip 091117
By=0.714;Cy=1.4;Ey=-0.2;
Br=0.852;Cr=2.3;Er=-2.75;

% Bx=.0822;Cx=1.65;Ex=-10;
% By=.348;Cy=1.799;Ey=-.184;
% By=0.714;Cy=1.4;Ey=-0.2;  %as in fdlat
% Bx=.2;Ex=.9;
Dx=1;Dy=1;

sx=0:100;
Fx=pacejka(sx,0,Bx,Cx,Dx,Ex,0);
aad=(0:30)/2;
Fy=pacejka(aad,0,By,Cy,Dy,Ey,0);
if nargin<1
    figure(1)
    subplot(211),plot(sx,Fx);ylabel('Fx')
    subplot(212),plot(aad,Fy);ylabel('Fy')
    % keyboard
end

%introduce Beckman combined slip model
issm=find(Fx==max(Fx));sxm=sx(issm);s=sx/sxm;
iaam=find(Fy==max(Fy));aadm=aad(iaam);a=aad/aadm;
Mx=[];My=[];
for ai=a
    Fxv=[];Fyv=[];
    for sj=s
        rho=sqrt(sj^2+ai^2);
        if rho==0
            Fxp=0;Fyp=0;
        else
            %                 om=Vx*(rho*SRm+1)/wheel.R;
            %                 Vy=Vx*tan(rho*aam);
            sxj=sxm*sj;aadi=ai*aadm;
            Fx=pacejka(rho*sxm,0,Bx,Cx,Dx,Ex,0);
            Fy=pacejka(rho*aadm,0,By,Cy,Dy,Ey,0);
            %                 F=ftyre(z,Vx,Vy,Vz,om,tyreparams,wheel,susp);
            Fxp=sj/rho*Fx;
            Fyp=ai/rho*Fy;
        end
        Fxv=[Fxv Fxp];
        Fyv=[Fyv Fyp];
    end
    Mx=[Mx; Fxv];
    My=[My; Fyv];
end
if nargin<1
    figure(2),surf(s*sxm,a*aadm,Mx),xlabel('ss (x100)'),ylabel('aa (deg)'),zlabel('Fx')
    figure(3),surf(s*sxm,a*aadm,My),xlabel('ss (x100)'),ylabel('aa (deg)'),zlabel('Fy')
end
return

function Y=pacejka(X,Sx,B,C,D,E,Sy)
% Pacejka magic formula
x=X+Sx;
Y=D*sin( C*atan( B*x - E*( B*x - atan(B*x)))) +Sy;
return
