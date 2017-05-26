function [loc,air,bm,act,wind,X]=db_droni(droni_xls,lat0,lon0)

%141121: airship model data for droni partially using worksheet
%
%original based on Aurora/dbm60
%
%----------------- STEPS ARE :
%-1-atmosphere and gas constants: air
%-2-experimental site local data: loc
%-3-blimp envelope geometric data: bm
%-4-mass and inertial data: bm
%-5-aerodynamic model data: bm
%-6-actuators: act
%-7-wind and turbulence model: wind
%-8-setup and simulation: X
%
%141125: still:
%> propulsion alternatives in forbm65 but one is chosen in xls
%> tail characteristics TBD

%0: assume droni configuration is to be retrieved from xls sheet
global Txt Num
if nargin<1,droni_xls='bm65_DRONI.xls';end
droni_sheet='DRONI';
% [status, sheetNames] = xlsfinfo(droni_xls)
[Num,Txt]=xlsread(droni_xls,droni_sheet);

modelst='preDRONI';
fprintf('\n Loading data for %s airship modelling using %s\n',modelst,droni_xls);
fprintf('(141121: adapted from dbm60/bm64//bm645)\n');

%1: air and constants -----------------------------------------
global g deg air
deg    =   pi/180;             % rad
g      =   9.80665;            % standard gravity acceleration (m/s2)

air=set_air;

%2: site specific ---------------------------------------------
if nargin<4
    lat0d=-3.38;lon0d=-64.72;%-> use defaults for Tefe'
end

loc.TG_C=getsingle('T_C');
loc.PG=getsingle('P_kPa')*1000;
loc.RH=getsingle('RH_%')/100;
hG0=(air.Po-loc.PG)/g/air.rro;

fprintf('Ground conditions set to: (hG=%dm) P=%.0f hPa, T=%.1f degC, RH=%.0f%%\n',...
    round(hG0),loc.PG/100,loc.TG_C,loc.RH*100);

%non-ISA atm. computation
loc.hG = air.To/air.dTdh*(1-(loc.PG/air.Po)^(1/(air.kk+1))); %pressure altitude above SL
loc.TG = air.To + (loc.TG_C - 15);
air.Ts = air.To - air.dTdh*loc.hG;
air.Td = loc.TG - air.Ts;                %difference to ISA temp. assumed as constant
rr_dry = air.rro*(air.Ts/air.To)^air.kk*(air.Ts/loc.TG);

%140402 include humidity input
loc.rrG=getdensity(loc.PG,loc.TG,loc.RH);

% local corrections on gravity: jra 050804
loc.lat0d = lat0d;    loc.lon0d = lon0d; 

%compute local gravity
g0   =9.7803184;  A_g =0.0053024; B_g =0.0000059; C_g =-3.086e-6;
g    =g0*(1+A_g*sin(loc.lat0d*deg)^2-B_g*sin(2*loc.lat0d*deg)^2)+C_g*loc.hG;
clear g0 A_g B_g C_g
%or use WGS84 model: for Campinas: g=9.7862

%--------------------------------------------------------------------------
%-3-blimp envelope geometric data: bm

bm.L=getsingle('length');
bm.d=getsingle('diameter');
bm.af=getsingle('front semi-axis');

bm.at     = bm.L-bm.af;          % from max diameter to tail end 
bm.Vol    = 1/6*pi*bm.d^2*bm.L;  % envelope volume
bm.xB     = 3/8*bm.at+5/8*bm.af; %=4.90m CB location from nose
bm.b      = bm.Vol^(1/3);        % equivalent blimp scale length

%ballonet
ceiling=getsingle('h_max');
overPres=getsingle('overPres');
bm.rbmax=getsingle('ballonet radius');
bm.Vbmax  = (4/3*pi*bm.rbmax^3); % max. ballonet volume (spherical)
bm.mb=getsingle('ballonet env.'); %141122 added ballonet mass
xb_from_nose=getsingle('ballonet from nose');
bm.OCairb = [bm.xB-xb_from_nose 0 bm.d/2];     
bm.GGmin  = 1-bm.Vbmax/bm.Vol;   % He minimum volume ratio
%validation
Vol_xls=getsingle('Volume Vol'); if abs(Vol_xls-bm.Vol)>1e-4*bm.Vol,disp 'ERR: Vol';end
xB_xls=getsingle('CV from nose'); if abs(xB_xls-bm.xB)>1e-4*bm.xB,disp 'ERR: xB';end
Vbmax_xls=getsingle('Vol_b_G'); if abs(Vbmax_xls-bm.Vbmax)>1e-4*bm.Vbmax,disp 'ERR: Vbmax';end

if bm.rbmax==0
    fprintf('--> NO BALLONET is considered: overPressure %.2f kPa at %.0fm\n',overPres,ceiling);
else
    fprintf('with %.0f%% BALLONET: design ceiling %.0fm\n',(1-bm.GGmin)*100,ceiling);
end

%computation of jB = inertia matrix ratio of displaced air (ellipsoid) 
% final inertia is multiplied by buoyancy mass mB (JB=jB*mB)
d2     = bm.d/2;
jxB    = 2/5*d2^2;
jyB0   = 19/320*(bm.at^2+bm.af^2)+1/5*d2^2+13/160*bm.af*bm.at;
bm.jB     = diag([jxB jyB0 jyB0]);

%%%%%%%%% virtual inertia ratios % from Ely_990430 /upmass_.m
%virtual terms: mv=diag([mmx mmy mmz])*mB/Jv=jv*mB
mmx=.0865;  mmy=.865;      mmz=mmy;
mmix=0;     mmiy=0.61;     mmiz=mmiy;     mmixz=0;
ri2=(bm.L^2+bm.d^2)/20;
%virtual inertia ratios (Mv=mv*mB,Jv=jv*mB...)
bm.mv=diag([mmx mmy mmz]);
bm.jv=[mmix 0 -mmixz;0 mmiy 0;-mmixz 0 mmiz]*ri2;
bm.jpqr=zeros(3); bm.juvw=zeros(3);

%adapt for bm5x/dv0x compatibility 050616
% Ji= [Jxi 0 -Jxzi;0 Jyi 0;-Jxzi 0 Jzi];bm.Jf=Ji-bm.jB*rr*air.MHe/air.Mair*bm.Vol;
bm.j_bo=zeros(3); bm.Vbo=0; bm.OCbo=zeros(1,3);  %dv0x compatibility: outer bal.
bm.DPHe=300;                               %He valve maximum overpressure ????
bm.eey=0;bm.eez=00;                           %ballonet linear volume for AS800

rr=getsingle('rho wet air');
bm.hg=getsingle('gondola height');

%-4-mass and inertial data: bm

%141125:for now weighting mass as is given by xls, not forced
[bm.mf,bm.mHe,bm.OCf,bm.Jf,mw]=cgcalc(rr,air,bm);

%-5-aerodynamic model data: bm
bm=set_aero(bm);

%ZZZZZZZ still missing for DRONI ----------------
% bm.lt = 4.01 ;               % tail fins moment arm
bm.lt = getsingle('fins aero')-bm.xB ;               % tail fins moment arm

%-6-actuators: act
%actuators from bm45/bm46
act.vect_min  =-30*deg;              %minimum vectoring angle (down)
act.vect_max  =120*deg;              %maximum vectoring angle (up)
act.vect_rate = 3;                   %maximum vectoring rate (rad/s);
act.def_max   = 25*deg;              %deflections saturation (symmetrical)
act.def_rate  = 1.5;                 %deflection maximum rate (rad/s)
act.servodelay= 0.05;                %engines servo delay (s)

%retrieve prop location in x from nose
prop_FR=getarray(4,'motor+prop FR');
prop_RR=getarray(4,'motor+prop RR');
prop_RL=getarray(4,'motor+prop RL');
prop_FL=getarray(4,'motor+prop FL');
prop_loc=([prop_FR;prop_RR;prop_RL;prop_FL]*[zeros(1,3);eye(3)])';
prop_loc(1,:)=prop_loc(1,:)-bm.xB; %into ABC: translation
act.OCTi=diag([-1 -1 1])*prop_loc; %into ABC: rotation

%141230 motor Temp. model: according to U8-10 in ecalc.ch
act.Kth=[1 1.4 2 3.5]; % [poor, medium, good, excellent] cooling

%conversion from de,da,dr to d1,d2,d3,d4 jra 001212
fprintf('for now X tail is considered in simulink as in AS800/DIVA\n')
bm.dear2d4=[1 -1 -1 1;1 1 1 1;1 1 -1 -1]'; %060303 normalize X deflections

%empirical model is (fbm2x): TBC
%  dwtail=1/4/(1+(dv-ddv0)^2/ddv2);
%  dwtail_f=(1+sqrt(1+((TR+TL)>0)*(TR+TL)*dwtail/rr/Sprop/(Vt+.001)^2))^2/4;
act.ddv0 = -.7/3; act.ddv2 = act.ddv0^2;     %angular alignment and width

%-7-wind and turbulence model: wind
%added to choose aerod.model
X.aeromodel=0;%bm45
% X.aeromodel=1;%AS800WT still to correct 050725
% X.aeromodel=2;%AS800WT coefs for forces only (higher drag and side force)
% X.aeromodel=3;%051207: AS800W2-->bm55
% X.aeromodel=4;%060123: AS800W3-->bm55 ; based on AS800WT-JRA060121.xls

%added for development to enable direct force inputs
X.forceinput=0; %default is real actuator inputs

%added for turbulence new model bm551 060919
X.drydenmodel=0; %bm45 model
% X.drydenmodel=1; %bm554 MILF dryden model
% X.drydenmodel=2; %bm554 MILH dryden model
wind.Tsw=.1; %110628: sampling for Band.Lin. white noise

%added to correct thom error in kinematics
X.thom=2;

%-8-setup and simulation: X
%sampling rate for output to workspace (tout,inout,nlout)
X.Toutput=.1;

%simulation results are in follwing output matrices: tout, inout, nlout
% sampled at 1/Toutput=10Hz into workspace
% .tout is time vector
% .inout is the dynamic function inputs 
% .nlout is the dynamic function outputs

                                       % WIND inputs SETTINGS
wind.ws=0;     wind.wh=0;                        % initial constant wind (m/s) and heading (rad)
wind.dws=0;    wind.dwh=0;     wind.tws=30;           % for wind step at time tws (s) (polar additive)
wind.turb=0;   wind.sigg=3;                      % gust (turb=1 is yes), sigg is intensity (m/s)

X.mis=[-100 -100;100 100];

%%%%%%%% variables are setup for fbm.m
fprintf('Data computed (%s).\n',modelst);
disp '------------------------------------------'

return

%-----------------------
function [mf,mHe,OCf,Jf,mw,mair_G,GG]=cgcalc(rr,air,bm)

%-->> fix mass is everything except: He + ballonet(air+env)

%assumed as in first column
[m_P,row_P]=getsingle('propulsion');
[m_A,row_A]=getsingle('airship');
[m_R,row_R]=getsingle('robotics');
[mf_xls,row_f]=getsingle('fix mass');

%retrieve [m x y z] 
% where coord. are x-from nose, y-left, z-down from env.axis
M_env=getarray(4,'envelope');
[M_P,T_P]=getmatrix(row_P(1)+1,row_A(1)-1,4);
[M_A,T_A]=getmatrix(row_A(1)+1,row_R(1)-1,4);
[M_R,T_R]=getmatrix(row_R(1)+1,row_f(1)-1,4);

Mass_fix=[M_env;M_P;M_A;M_R];
mf=sum(Mass_fix(:,1));

OCf_from_nose=Mass_fix(:,1)'*Mass_fix(:,2:4)/mf; 
OCf = (OCf_from_nose - [bm.xB 0 0])*[-1 0 0;0 -1 0;0 0 1]; %in ABC

M_He_G=getarray(4,'helium');
M_bair_G=getarray(4,'ballonet air');
M_benv_G=getarray(4,'ballonet env.');

%141125:for now weighting mass as is given by xls, not forced
Mt_G=[Mass_fix;M_He_G;M_bair_G;M_benv_G];
mt=sum(Mt_G(:,1));
B=rr*bm.Vol;
mw=mt-B;
%----------------------------------------------------------------
mHe=M_He_G(1);
mair_G=M_bair_G(1);
GG=1-bm.Vbmax/bm.Vol;
mHe   = GG*rr*bm.Vol*air.MHe/air.Mair;                        % He mass 
Vair  = bm.Vol*(1-GG);                                % ballonet air volume
mair  = Vair*rr;                                   % ballonet air mass
if bm.Vbmax==0, rb=0;
else
    rb=bm.rbmax*Vair/bm.Vbmax;                %>>050107: assumed linear relation
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Computing OC - c.g (less air)
OCi_from_nose=Mt_G(:,1)'*Mt_G(:,2:4)/mt; 
OCi = (OCi_from_nose - [bm.xB 0 0])*[-1 0 0;0 -1 0;0 0 1]; %in ABC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
M=Mt_G(:,1);

x=bm.xB-Mt_G(:,2); y=-Mt_G(:,3); z=Mt_G(:,4);

Jx=(M.*(y.*y+z.*z)); 
Jy=(M.*(x.*x+z.*z)); 
Jz=(M.*(y.*y+x.*x)); 
Jxz=(M.*x.*z); 

% inertia principals 
data_He=getarray(3,'rho_He_ISA');
rr_He=data_He(3);

%141125: probably to include into xls also!!
Jx_g=.07;Jy_g=.15;Jz_g=.16;                 % gondola (2) < bm46
%---------

% --- Helium
% d2=d/2;
% VolB =2/3*d2^2*pi*L
% jxB =2/5*d2^2
% jyB0=19/320*(at^2+af^2)+1/5*d2^2+13/160*af*at
d2=bm.d/2;
jxB=2/5*d2^2;
jyB=19/320*(bm.at^2+bm.af^2)+1/5*d2^2+13/160*bm.af*bm.at;
Jx_He=jxB*rr_He*bm.Vol;
Jy_He=jyB*rr_He*bm.Vol;Jz_He=Jx_He;Jxz_He=0;

% --- envelope: symbolic, from sydv00
ss=getsingle('env.specific')/1000; %env.specific mass in kg/m2
Jx_e =3/8*pi^2*d2^3*bm.L*ss;
%130729: some sign error for Jy_e Jz_e???
Jy_e=-(-1/2*d2*bm.L*(bm.at-bm.af)^2*pi+1/128*d2*bm.L*(84*bm.af*(bm.af-bm.at)+25*bm.L^2+24*d2^2)*pi^2*ss);
Jz_e=Jy_e;Jxz_e=0;
%--------
        
Jxf=sum(Jx)+Jx_e+Jx_g;
Jyf=sum(Jy)+Jy_e+Jy_g;
Jzf=sum(Jz)+Jz_e+Jz_g;
Jxzf=sum(Jxz);   

Jf=[Jxf 0 -Jxzf;0 Jyf 0;-Jxzf 0 Jzf]; 

return

%-----------------------
function bm=set_aero(bm)
global deg
%%%%%%%%% aerodynamic coefficient derivatives <bm45
% possibly to define in aerodynamic model function

% (some of them could be related to geometrical formulae: McCormick)
%from Ely_990430 /aerodd_.m (all C in /rad)
%revisited on 001214 jra/xcoef(ely_3049)
bm.Cd0=0.0326;      bm.Cdi=1.25;
bm.Cl0=-.0023;      bm.Cla=.016/deg;    bm.CYb=-bm.Cla;  %in CYb sign changed 001214
bm.CM0=.0017;       bm.CMa=.0057/deg;   bm.CNb=-.01/deg;
bm.Clde=.005/deg;   bm.CYdr=bm.Clde;
bm.CMde=-.0045/deg; bm.CNdr=-.0054/deg; bm.CLda=-.006/deg;
bm.CMab=-0.458;     bm.CMba=-0.3;                  %2nd.added by jra 001214
bm.CLb =  4*0.1500; bm.CMb =   -0.1200;            %added by ely 991025/revisited jra001224

% damping dimensions corrected 990823/990923: cxx=Cxx*L/2/V (s/rad)
bm.cmq=-0.256;      bm.cnr=-.258;       bm.clp=-.012;

% due to missing experimental data, tail surfaces coef. model based on MC
deduce_tail_aero=1;if deduce_tail_aero
    %as800b tail fins data 050621
    cfin=(.85+.4+.45+.34)/2;bfin=.67;cflap=(.4+.34)/2;xn2fin=7.9;

    %150128 deduce approximate tail arm from xls
    fins=getarray(4,'tail fins'); xn2fincg=fins(2);
    xn2fin=xn2fincg-cfin/2;
    %150128-> how to correct virtual mass terms ? TBC
    %---
    
    Sfin=cfin*bfin;
    clat=.1/deg;nn=.7;
    ttf=acos(2*cflap/cfin-1);taut=1-(ttf-sin(ttf))/pi;
    S=bm.Vol^(2/3);
    bm.Cldd=Sfin/S*clat*taut*nn; %local lift coef. per fin
    nn_x=.75;%nn_x=1;
    bm.Clde=2*sqrt(2)*bm.Cldd*nn_x;
    lf=xn2fin+cfin/4-bm.xB;  %should be equal to lt ??? 050724
    bm.CMde=-lf/bm.b*bm.Clde;
    bm.CYdr=bm.Clde;
    bm.CNdr=1.2*bm.CMde;bm.CLda=4/3*bm.CMde; %assume similarity w/YEZ
    fprintf('Tail surf. coefs computed: Clde=%.4f, CMde=%.4f (/rad)\n',bm.Clde,bm.CMde);
end
return

%-----------------------------
function air=set_air
%definitions of general constants
air.RR     =   8.314;              % perfect gas constant (SI)
air.MHe    =   0.004;              % mass of a mole of He (kg)  
air.Mair   =   0.028964;           % mass of a mole of air (kg)
air.dTdh   =   0.0065;             % ISA temperature gradient (K/m)
air.To     = 288.15;               % ISA SL temperature (K) = 15.C
air.Po     =  101325;              % ISA SL pressure (Pa)
air.rro    =   1.225;              % ISA SL density (kg/m3)
air.kk     =   4.2558797;          % ISA temp to density exponent

%140402 NEW add humidity input < airdensity-140402.html
air.Rair=287.058; air.Rvap=461.495; %perfect gas constants
return



% site specific ---------------------------------------------


function rho_wet=getdensity(P_s,T_K,RH,buck)

% function rho_wet=getdensity(P_s,T_K,RH)
%
%140402: compute wet air density including humidity input
%
% P_s   static pressure in Pa
% T_K   temperature in K
% RH    relative humidity 0 < RH < 1
% <buck> if defined then use Buck's formula
%
% rho_wet : final density (kg/m3): rho_wet = rho_dry + rho_vap
%
%< airdensity-140402.html or http://wahiduddin.net/calc/density_altitude.htm
%< alt. andreasHB/Buck's formula ?ref
%
%--> water vapor is lighter than air (18g vs 29g g/mol): density is reduced

global air
if length(air)==0; fprintf('ERR: first air must defined\n');end

if nargin==3 %< airdensity-140402.html
    Psat=601.78*10^((7.5*T_K-2048.625)/(T_K-35.85));
else %< alt. andreasHB/Buck's formula
    Psat=611.21*(1.0007+3.46e-8*P_s)*exp(17.502*(T_K-air.To)/(240.97+T_K-air.To));
end
Pv=RH*Psat;
rho_dry=(P_s-Pv)/air.Rair/T_K;
rho_vap=Pv/air.Rvap/T_K;
rho_wet=rho_dry+rho_vap;
return

%----------------------- retrieve data from xls -----------------
function [var,found]=getsingle(name,var)
%retrieve single parameters from xls
global Txt Num
[m,n]=size(Txt);
found=0; for r=1:m
    for c=1:n,
        if strfind(Txt{r,c},name),found=[r,c];break,end
    end
end
if found, value=Num(found(1),found(2)); else value=[]; end
if ~isnan(value), var=value;
else if nargout==1,fprintf('ERR: missing "%s"\n',name);end
    if nargin<2,var=[];end
end 
%----------
function [var,found]=getarray(len,name,var)
%retrieve array parameters from xls
global Txt Num
[m,n]=size(Txt);
found=0; for r=1:m
    for c=1:n,
        if strfind(Txt{r,c},name),found=[r,c];break,end
    end
end
if found,if found(2)+len-1>size(Num,2),found=0;end,end
if found, value=Num(found(1),found(2)+(0:len-1)); else value=[]; end
if ~isnan(value), var=value;
else if nargout==1,fprintf('ERR: missing %s\n',name);end
    if nargin<2,var=[];end
end 
%----------
function [var,name]=getloc(len,found)
%retrieve parameters from xls knowing loc
global Txt Num
if found(2)+len-1>size(Num,2),var=[];name=[];return
else 
    name=Txt{found(1),found(2)};
    var=Num(found(1),found(2)+(0:len-1));
end
%----------
function [M,T]=getmatrix(top,low,len)
M=[];T=[];for r=top:low
    [value,item]=getloc(4,[r,1]); 
    if length(item)>0
        M=[M; value];T=strvcat(T, item); end
end
return
