
% October 2008 
% Write your name here ! 
% Plot of Characteristics for Dugoff tire model 
% -------------------------------------------------------- 
clear all 
clc 
% --------------------------------------------------------- 
% Enter tire data for Dugoff tire model 
% --------------------------------------------------------- 
% Enter velocity component in the wheel plane 
%      Ukph=?????; % [kph] 
      U=15; % [m/s] 
% Enter longitudinal slip 
     S=0.0; 
% Enter longitudinal tire stiffness [N/unit slip] 
Bx=.0822;Cx=1.65;
Cx=Bx*Cx;
Cx=Cx*100;
% Enter friction reduction factor [s/m] 
     As=0.0115; 
% Enter peak static tire/road friction coefficient 
     mu0=1; 
   
% Enter tire cornering stiffness - in [N/rad] 
     By=0.714;Cy=1.4;
     Ca=By*Cy;
     Ca=Ca*180/pi;
% --------------------------------------------------      
   i=0;  Fx=[];  SS=[];
   alfa=1e-7; % [rad] 
    for S=0:0.005:1;
    i=i+1;
    SS(i)=S;
        z=(1+S)/(2*sqrt(((Cx*S)^2)+(Ca*tan(alfa))^2));
        if z<1 
            f=z*(2-z); 
        else 
            f=1; 
        end 
        Fx(i)=Cx*S*f/(1+S); 
    end
     figure(5)
     subplot(211),plot(SS*100,Fx);ylabel('Fx')
% ---------
      i=0;  Fa=[];  SS=[];
      S=0;
     for alfa=0:0.05:15;
         alfa=alfa*pi/180; % [rad] 
    i=i+1;
    aa(i)=alfa*180/pi;
        z=(1+S)/(2*sqrt(((Cx*S)^2)+(Ca*tan(alfa))^2));
        if z<1 
            f=z*(2-z); 
        else 
            f=1; 
        end 
        Fa(i)=Ca*tan(alfa)*f/(1+S); 
    end
     figure(5)
     subplot(212),plot(aa,Fa);ylabel('Fa')
 % ---------    
     
Mx=[];My=[];
s=0:0.005:1; s=s*100
a=0.0001:0.5:15; 

for ai=0.0001:0.5:15;
    ai=ai*pi/180; % [rad] 
    Fxv=[];Fyv=[];
    for sj=0:0.005:1;
     
        z=(1+sj)/(2*sqrt(((Cx*sj)^2)+(Ca*tan(ai))^2));
        if z<1 
            f=z*(2-z); 
        else 
            f=1; 
        end 
        Fyp=Ca*tan(ai)*f/(1+sj); 
        Fxp=Cx*sj*f/(1+sj);    
                     
        Fxv=[Fxv Fxp];
        Fyv=[Fyv Fyp];
    end
    Mx=[Mx; Fxv];
    My=[My; Fyv];
end
% if nargin<1
    figure(6),surf(s,a,Mx),xlabel('ss (x100)'),ylabel('aa (deg)'),zlabel('Fx')
    figure(7),surf(s,a,My),xlabel('ss (x100)'),ylabel('aa (deg)'),zlabel('Fy')
% end
 
 
 
 
 
 
 
j=0; % Counter for normal tire load in [N] 
angmax=15
for Fz=2000:2000:6000 
    i=0; % Counter for slip angle 
     j=j+1; 
    for alfa=0.0001:0.05:angmax;  % [deg] 
        alfa=alfa*pi/180; % [rad] 
        i=i+1; 
        %Vs=U*sqrt(S^2+(tan(alfa))^2); 
        %mu=m0*(1-As*Vs); 
        z=(1+S)/(2*sqrt(((Cx*S)^2)+(Ca*tan(alfa))^2))
        if z<1 
            f=z*(2-z); 
        else 
            f=1; 
        end 
        %   Fx=Cl*S*f/(1-S); 
        Fyf(i,j)=Ca*tan(alfa)*f/(1+S); 
    end 
end 
% -------------------------------------------------- 
% Plot Characteristics 

figure(4)
alfa=0.0001:0.05:angmax; 
plot(alfa,Fyf(:,1),alfa,Fyf(:,2),alfa,Fyf(:,3),'linewidth',2) 
grid on 
title('Dugoff Tire Cornering Force Characteristics','fontsize',12,'fontweight','b') 
legend('F_z = 2 kN','F_z = 4 kN','F_z = 6 kN','Location','NorthWest') 
xlabel('Slip angle, \alpha(^o)','fontsize',12,'fontweight','b') 
ylabel('F_y/F_z','fontsize',12,'fontweight','b') 
 