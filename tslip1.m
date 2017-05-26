%tlat.m/tslip.m
%090924/091117
%
% script to model lateral sliding tyre dynamics
%091117: added longitudinal part for combined model
%
%--> setup script for sdslip/fdslip
%091118: added 
% -kinematics in sdslipxy
% -trim+linearization
% -guidance control assuming no sliding in sdslipxy
% -lon regulator at trim
%
%still:
% -may add low level close loop to regulate sliding and track yaw rate
% -may add lon. tracking control
% -add slip/slide estimator ? added 091124

% vehicle  parameters
global acar b c d m Jr g mu deg cont
cont=1;

%a=1.5;b=2;c=.9;d=.5;
acar=.7;b=.6;c=.59;d=0.45;
% d e' a altura do c.g. acima do solo, utilizado para deduzir a forca centrifuga
%a=.5;b=.8;

%m=1200;
m=800;

Jr=1/2*m*1.5^2;
g=9.81;deg=pi/180;
% tyres model in fdslip as in ttyre2 < Milliken / Brian Beckman
% (model parameters should be defined in setup script)
mu=1.7/2;  %reduced friction coef. to lower unstability
% mu=1.7/10;  %further reduced friction coef. to lower unstability
mu=0.8;

%added for tslip/fdslip
global Jw Rw bx bom
Rw=0.25;Jw=2.5*10*.2^2; %wheel spinning inertia
bx=1;bom=.003;

% tyres model as in ttyre_Pacejka < Milliken / Brian Beckman
% global sxm aadm
ttyre_Pacejka(1); %same model curves as in fdslip
% disp(sprintf('sxm=%g, aadm=%g deg',sxm,aadm))

% add longitudinal to vertical forces coupling
global Fx0
Fx0=0;
% syms Fz12 Fz34 Fx a b d mg real
% Emom=a*Fz12-b*Fz34-d*Fx;
% Efor=Fz12+Fz34-mg;
% s0=solve(Emom,Efor,'Fz12','Fz34');
% [s0.Fz12;s0.Fz34]

global traccoef % 0:back/1:front;
traccoef=0;

%--------- prepare data for simulation
tsim=150;

ui=.01;Vi=[ui 0 0 [1 1 1 1]*ui/Rw];
xi=0;yi=0;psii=0;
% sdslip
% return

%--------- trim and linearize dynamics at selected speed u0
% XU=[ u v r [om] dd T]
u0=15; %limit of stability for mu=0.17 (reduced case)
% u0=10;
% u0=5; %reduced slip
u0=3.5;   % u0=4.5 com mu=0.4
        % u0=6.3 com mu=0.8 (SMC OK, mas PD NAO, LQR NAO)
        % u0=1.5 com mu=0.05 (SMC OK, PD OK, LQR NAO, apenas com mu>0.1)
global uref
uref=u0;
%initialize trimming
v0=0;r0=0;om0=u0/Rw;dd0=0;T0=44*(u0/10)^2;
XU0=[u0 v0 r0 om0*[1 1 1 1] dd0 T0]';
xxX=eye(7,15+4);
dXdt0=xxX*fdslip12(XU0);
XU00=XU0;

%find trim condition
dd=.000001;
for rep=1:5
    ff=[];kkv=[4:7 8 9];
    for kk=kkv
        XU1=XU0+((1:9)'==kk)*dd;
        dXdt1=xxX*fdslip12(XU1);
        ff=[ff (dXdt1-dXdt0)/dd];
    end
    % f1=f0+Df*(x1-x0)=0  --> x1=x0-Df^-1*f0
    XU0(kkv)=XU0(kkv)-pinv(ff)*dXdt0
    dXdt0=xxX*fdslip12(XU0)
    if max(abs(dXdt0))<1e-4;break;end
end
disp 'final progress'
dXdt0'
if max(abs(dXdt0))>.001;
    disp(sprintf('!!! did not succeed to reach trim condition for u=%gm/s',u0));
    disp '-------LEAVING------'
    return
else
    disp 'trim condition XU0=[u v r om1 om2 om3 om4 dd T]'
    XU0'
end

%linearize dynamics
ff=[];kkv=1:9;
for kk=kkv
    XU1=XU0+((1:9)'==kk)*dd;
    dXdt1=xxX*fdslip12(XU1);
    ff=[ff (dXdt1-dXdt0)/dd];
end
ff

s7=ss(ff(1:7,1:7),ff(1:7,8:9),eye(7),zeros(7,2))
s3=modred(s7,4:7)
% @ 25m/s -----------------------------------
% a = 
%                 x1           x2           x3
%    x1     -0.07361  -2.226e-006  -8.498e-006
%    x2   1.756e-020         -1.4          -25
%    x3  -2.774e-017        -1.12       -3.672
%  
% b = 
%                 u1           u2
%    x1  -7.538e-005     0.003243
%    x2        19.29  -3.338e-021
%    x3        42.67   2.168e-019
damp(s3)
%  Eigenvalue      Damping     Freq. (rad/s)  
%  -7.36e-002     1.00e+000      7.36e-002    
%   2.88e+000    -1.00e+000      2.88e+000    
%  -7.95e+000     1.00e+000      7.95e+000    

 if 0
    X=[];%ABM=[];k=1;
    for u=1:1:6
        v=0;r=0;dd=0*deg;
        dvrdt=fdlat([u,v,r,dd]);
        deps=.0001;
        sprintf('linearized dynamics at u=%d m/s X=[v r] U=[dd]',round(u))
        AB=[fdlat([u,v+deps,r,dd]) fdlat([u,v,r+deps,dd]) fdlat([u,v,r,dd+deps])]/deps
        X=[X [u;eig(AB(1:2,1:2))]];%ABM(:,:,k)=AB;k=k+1;
    end
    X
    % @ 25m/s -----------------------------------
    damp(AB(:,1:2))
    %  Eigenvalue      Damping     Freq. (rad/s)  
    %   2.85e+000    -1.00e+000      2.85e+000    
    %  -7.99e+000     1.00e+000      7.99e+000    
end

[a3,b3]=ssdata(s3);
Ts=.01*5;
k3=lqrd(a3,b3,eye(3),diag([1 .01/1000]),Ts)

%set simulation inputs
Vi=XU0(1:7);dd0=0;T0=XU0(9);

showcl=0;showol=1;showmis=2;showcle=3;
usecase=showmis;
% usecase=showol;
%-------- here may define a low level close loop to regulate sliding
% usecase=showcl;
% usecase=showcle;

switch usecase
    case showcl
        a2=a3(2:3,2:3);b2=b3(2:3,1);
        k2=lqrd(a2,b2,eye(2),1,Ts);
        c2=[0 1];
        dcg2=-c2*inv(a2-b2*k2)*b2; f2=inv(dcg2)
        tstepdd=5;
        dstepdd=.1*deg;
        R=50;dstepr=u0/R;
        tsim=20;
        sdslip1r
        return
    case showcle        
        a2=a3(2:3,2:3);b2=b3(2:3,1);
        k2=lqrd(a2,b2,eye(2),1,Ts);
        c2=[0 1];
        dcg2=-c2*inv(a2-b2*k2)*b2; f2=inv(dcg2)
        
        % === side slip estimator =============== 091124
        %
        % how do you estimate slipping (for wheels?), what sensors ?
        %
        % dxdt=ku*T+w
        % dydt=u*syy+v=u*sin(yy+aa)
        % dyydt=kr*u*dd+n
        % dwdt=ew
        % daadt=eaa
        % dndt=en
        % X=[x y yy w aa n] U=[T dd] W=[ew eaa en]
        ku=0.22; %for u0=10
        kr=0.33; %for
        AX=[0 0 0 1 0 0; 0 0 u0 0 u0 0; 0 0 0 0 0 1; zeros(3,6)];
        BX=[ku 0; 0 0; 0 kr*u0; 0 0; 0 0; 0 0];
        GX=[zeros(3); eye(3)];
        % Y=[x yy lam]
        CX=[eye(1,6);0 0 1 0 0 0; 0 0 1 0 1 0];DX=zeros(3,2);
        Qn=1e-2*eye(3);Rn=1e-4*eye(3);
%         Qn=1e-1*eye(3);Rn=1e-4*eye(3);
        % remove y state
        AX(2,:)=[];BX(2,:)=[];GX(2,:)=[];AX(:,2)=[];CX(:,2)=[];
        Kest=kalman(ss(AX,[BX GX],CX,[DX zeros(3)]),Qn,Rn);
        [ae,be,ce,de]=ssdata(Kest);
        
        tstepdd=5; tstepdd1=15; dstepdd=.1*deg;
        R=50;dstepr=u0/R;
        tsim=20;
        gointeractive=0;if gointeractive
            sdslip1re
            return
        else
            y1=sim('sdslip1re');
            % at the end of simulation, run following to plot results....
            %------------
            % yout=[t, u v r oms, d[u v r om], F sig aa, dd T,...
            %     xe yye we aae nre, x y yy] (1+7+7+3*7+2+5+3)
            xy=yout(:,1+7+7+3*4+2+5+(1:2));
            figure(1),plot(xy(:,2),xy(:,1)),title('xy path [m x m]')
            tout=yout(:,1);dt=(tout(2)-tout(1));
            vxy=[[0 0];diff(xy)]/dt;
            v=sqrt(sum(vxy'.^2))';v(1)=mean(v(2:10));
            figure(2),subplot(211),plot(tout,vxy),ylabel('v_x,v_y [m/s]')
            subplot(212),plot(tout,v-u0),ylabel('|v_{err}| [m/s]')
            aas=yout(:,1+7+7+2*4+(1:4)); % model wheels side slip angles
            aae=yout(:,1+7+7+3*4+2+3+1); %estimated vehicle side slip
            lam=atan2(vxy(:,2),vxy(:,1)); %path angle
            yy=yout(:,end); % model yaw angle
            figure(3),subplot(211),plot(tout,[-aas aae]),title('(-1*) slip angles & est.aa')
            subplot(212),plot(tout,[aae rem(lam-yy,pi)]),legend('est.aa','path slip angle')
            % plot xy plus yaw
            Vpts=[a -c;a c;-b c;-b -c;a -c];
            figure(10),clf,plot(xy(:,2),xy(:,1)),hold on
            ddt=10;for kk=1:ddt:length(tout)
                R=[cos(yy(kk)) sin(yy(kk));-sin(yy(kk)) cos(yy(kk))];
                ptsk=5*Vpts*R;ptsk(:,1)=ptsk(:,1)+xy(kk,1);ptsk(:,2)=ptsk(:,2)+xy(kk,2);
                plot(ptsk(:,2),ptsk(:,1),'g')
            end
            hold off
            mosaic(3)
            return
        end
    case showol
        
        tstepdd=5;
        dstepdd=.001*deg;
        tsim=20;
        sdslip
        %         return
        % plot curvature gain: 1/R/dd vs time
        fig=1;for u0=[2 5 8 11 14 17 20]
            T0=44*(u0/10)^2;Vi(1)=u0;
            sim('sdslip1');
            figure(fig),fig=fig+1;plot(yout(:,1),yout(:,4)./yout(:,2)/dstepdd)
        end
        mosaic(fig-1)
        return
    case showmis
        
        % show mission tracking guidance w/o SC
        
        % initial variables and mission setup for sdlatxy
        ddstep=5*deg;
        % constant lon.speed
        u=u0; %: sideslip: unstable but tracks!
        
        %guidance controler assuming no sideslip
        kg=.025/ddstep/5;
        kg=0.76;
        
        Ag=[0 u;0 0];Bg=[0; kg*u];
        Kgcar=lqr(Ag,Bg,1.0*eye(2),10);
        xi=0;yi=1;psii=0.0;
        
        mu_0=0.8;
        Ag2=[-1.6/m/u (b*mu_0-acar*mu_0/m/u)-u; (b*mu_0-acar*mu_0/Jr/u) -(acar*acar*mu_0+b*b*mu_0/Jr/u)];
        Bg2=[mu_0/1 ; acar*mu_0/1];
        
        Ag1=[  -20.59          -4;
                -10.8      -7.686];
            
        Bg1=[38; 43.6];    
        Kgcar=lqr(Ag1,Bg1,[10 0;0 10],10);
        
        global mis imis dmis circ
        
        test45deg=0;if test45deg
            % straight segments
            circ=0;
            mis=[0 0;100 0;140 50;140 100;100 140;0 140;0 0];
            tsim=250/u;
        else
            % circle segments with smooth path
            circ=1;
            %mis=[0 0;100 0;100 100;0 100;0 0];
            
%             raio=15;
%             dmis=raio/5;
            mis=[0 0 0;
                10 0 pi/2;
                15 5 pi/2;
                10 10 0;
                0 10 pi/2;
                -5 5 pi/2;
                0  0   0];
            
           % tsim=(50*1.2*pi+2*100)/u;
           % tsim=3*2*(raio*1.2*pi+2*2*raio)/u;
           
        end
        imis=1;%dmis=10;
        gointeractive=1;if gointeractive
            sdslip1xy
            return
        else
            y1=sim('sdslip1xy');
            ny=size(yout,2);tout=yout(:,1);nt=length(tout);
            xy=yout(:,ny-3+(1:2));
            yy=yout(:,end); % model yaw angle
            % plot xy plus yaw
            Vpts=[a -c;a c;-b c;-b -c;a -c];
            figure(10),clf,plot(xy(:,2),xy(:,1)),hold on
            ddt=10;
            for kk=1:ddt:length(tout)
                R=[cos(yy(kk)) sin(yy(kk));-sin(yy(kk)) cos(yy(kk))];
                ptsk=5*0.2*Vpts*R;
                
                ptsk=[0 0; 1 0];
                ptsk=ptsk*R;
                
                ptsk(:,1)=ptsk(:,1)+xy(kk,1);
                ptsk(:,2)=ptsk(:,2)+xy(kk,2);
                
                plot(ptsk(:,2),ptsk(:,1),'r')
                
                plot(ptsk(:,2),ptsk(:,1),'r')
            end
            hold off
            return
        end
end