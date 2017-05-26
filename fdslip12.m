function out=fdslip(inp)

% function duvromdt=fdslip(inp)
%
% dynamic function for slipping tyres car dynamics
% adapted version to include both lon. and lat. slipping
%
% input  =[ u v r [om] dd T]
% output =[ dXdt [Fx] [sig] [aa]]
% dXdt   =[ dudt dvdt drdt [domdt]]
%
% u : longitudinal forward speed
% v : side speed (to vehicle right)
% r : yaw rate
% om: wheels angular speed (4)
% dd: equivalent steering angle
%
% outputs: state derivatives, for all wheels lon.force/slip ratio/sideslip angle
%
%090924/091117
%100209: fdslip12: correction on vertical force due to

u=inp(1);v=inp(2);r=inp(3);om=inp(4:7);
% same steering angle assumed on front wheels (an ackerman model could be used)
dd=inp(8);
% distribution of torque on wheels is defined here
T=inp(9);
% Tv=[T T T T]'/4;  % equal for all wheels
% Tv=[0 0 T T]'/2;  % both back wheels
% Tv=[T T 0 0]'/2;  % both front wheels
global traccoef % 0:back/1:front; 100204
Tv=T/2*[[1 1]*traccoef [1 1]*(1-traccoef)]';

[FxyN,Fxv,sig,aa]=fdslip0(u,v,r,om,dd);

global m Jr Jw Rw bx bom

dudt=1/m*( FxyN(1) -bx*u*abs(u));
dvdt=FxyN(2)/m-u*r;
drdt=FxyN(3)/Jr;
domdt=1/Jw*( Tv -Fxv*Rw -bom*om.*abs(om));%if u>5;keyboard,end
duvromdt=[ dudt; dvdt; drdt; domdt];

out=[ duvromdt; Fxv; sig; aa];

return

function [FxyN,Fxv,sigv,aav]=fdslip0(u,v,r,om,dd)

% function to compute forces acting on tyres 
% knowing 2D velocity and wheels rotation speed, and steering angle
%
%Pacejka Fx, Fy and Mz (N) curves from Milliken

global acar b c d m g mu Rw cont

%dd1=dd;dd2=dd;
% Ackerman model for front wheels (dd1,dd2)

     L=acar+b;
     dd1 = atan(L*tan(dd)/(L+c*tan(dd)));
     dd2 = atan(tan(dd)*L/(L-c*tan(dd)));

if abs(u)<.01;aa=zeros(1,4);else
    aa(1)=atan((v+acar*r)/(u+c*r))-dd1; %FL
    aa(2)=atan((v+acar*r)/(u-c*r))-dd2; %FR
    aa(3)=atan((v-b*r)/(u+c*r)); %BL
    aa(4)=atan((v-b*r)/(u-c*r)); %BR
end

% wheels longitudinal speeds
uv(1)=(u+c*r)*cos(dd1);
uv(2)=(u-c*r)*cos(dd2);
uv(3)=u+c*r;
uv(4)=u-c*r;

% wheels vertical forces
Fz(1)=b/(acar+b)/2*m*g;
Fz(2)=b/(acar+b)/2*m*g;
Fz(3)=acar/(acar+b)/2*m*g;
Fz(4)=acar/(acar+b)/2*m*g;

%-----------------------------------------
% add longitudinal to vertical forces coupling CORR100209
global Fx0 %forward longitudinal force at cg
Fz(1:2)=Fz(1:2)+Fx0*d/(acar+b)/2;
Fz(3:4)=Fz(3:4)-Fx0*d/(acar+b)/2;
% add lateral to vertical forces coupling CORR100204/100209
Fy0=-m*u*r; %centrifugal force to right
Fz0=Fy0*d/2/c;
Fz(1)=Fz(1)-Fz0*b/(acar+b);
Fz(2)=Fz(2)+Fz0*b/(acar+b);
Fz(3)=Fz(3)-Fz0*acar/(acar+b);
Fz(4)=Fz(4)+Fz0*acar/(acar+b);
%cancel negative vertical force 100209
Fz=max([zeros(1,4);Fz(:)']);
%-----------------------------------------

% tyres model as in ttyre2 < Milliken / Brian Beckman
global Bx Cx Ex By Cy Ey Br Cr Er
% Bx=.0822;Cx=1.65;Ex=-10; %added for combined slip 091117
% By=0.714;Cy=1.4;Ey=-0.2;
% Br=0.852;Cr=2.3;Er=-2.75;
%for combined slip: maximum abscissae: sxm in %, aadm in deg
global sxm aadm
aam=aadm*pi/180;sm=sxm/100;

Fx=zeros(1,4);Fy=zeros(1,4);N=zeros(1,4);
Fx1=zeros(1,4);Fy1=zeros(1,4);N1=zeros(1,4);
sig=zeros(1,4);
% for i=1:4
%     aab=Cy*tan(aa(i))*15; %corrected side slip angle
%     V=sqrt(uv(i)^2*(1+aa(i)^2));if V<.001;V=.001;end
%     sig(i)=(om(i)*Rw-uv(i))/V; %slip ratio
%     a1=aab/aam;s=sig(i)/sm;rho=sqrt(a1^2+s^2);
%     if rho<1e-9;rho=1e-9;end
%     sx=rho*sxm;aab=rho*aam;
%     Dx=mu*Fz(i);
%     Fx1(i)=pacejka(sx,0,Bx,Cx,Dx,Ex,0);
%     Dy=mu*Fz(i);
%     Fy1(i)=pacejka(aab,0,By,Cy,Dy,Ey,0);
%     Dr=0.51*mu*Fz(i)*c;
%     N1(i)=pacejka(aab,0,Br,Cr,Dr,Er,0);
%     
%     Fx(i)=Fx1(i)*s/rho;
%     Fy(i)=Fy1(i)*a1/rho;N(i)=N1(i)*a1/rho;
% end
% keyboard

muu=mu; %+2*0.5*[0.2 -0.2 0.2 -0.2]

aab=Cy*tan(aa)*15; %corrected side slip angle
V=sqrt(uv.^2.*(1+aa.^2));V=max(V,.001*[1 1 1 1]);
sig=(om'*Rw-uv)./V; %slip ratio
a1=aab/aam;s=sig/sm;rho=sqrt(a1.^2+s.^2);
rho=max(rho,1e-9*[1 1 1 1]);
sx=rho*sxm;aab=rho*aam;
Dx=muu.*Fz;
Fx1=pacejka(sx,0,Bx,Cx,Dx,Ex,0);
Dy=muu.*Fz;
Fy1=pacejka(aab,0,By,Cy,Dy,Ey,0);
Dr=0.51*muu.*Fz*c;
N1=pacejka(aab,0,Br,Cr,Dr,Er,0);

Fx=Fx1.*s./rho;
Fy=Fy1.*a1./rho;N=N1.*a1./rho;
% figure(1),subplot(211),plot(aa/deg,Fy),subplot(212),plot(aa/deg,N)

Fxcg=sum(Fx.*cos([dd1 dd2 0 0]))-sum(Fy.*sin([dd1 dd2 0 0]));
Fycg=sum(Fy.*cos([dd1 dd2 0 0]))+sum(Fx.*sin([dd1 dd2 0 0]));
Ncg=sum(N)-b*(Fy(3)+Fy(4))+acar*cos(dd)*(Fy(1)+Fy(2))...
    +c*cos(dd)*(Fx(1)-Fx(2))+c*(Fx(3)-Fx(4));
FxyN=[Fxcg; -Fycg; -Ncg];
Fxv=Fx';
% if u>5;keyboard,end
sigv=sig';aav=aa';
return

function Y=pacejka(X,Sx,B,C,D,E,Sy)
% Pacejka magic formula
x=X+Sx;
Y=D.*sin( C*atan( B*x - E*( B*x - atan(B*x)))) +Sy;
return

