function fout=fbm65_kin(inp)

%function [output]=fbm65([x;U])
%
% (inp ==> [x,U], state, inputs)
%
% fout=[state other_outputs] (13+22)
% state=[u v w p q r PN PE PD e0 e1 e2 e3] (13)
% gust_state =[ug xvg1 xvg2 xwg1 xwg2 pg xqg xrg]  (8)
% other_outputs=[u v w p q r dN dE dD e0 e1 e2 e3 Vt [aa bb ff tt yy] ua va wa] (22)

%%%%%%%%%%%%%%%% variables and parameters
x=inp(1:21);                % (21) states
U=inp(22:32);               % (5+6) inputs (wind (5) and control (6))

u=x(1); v=x(2); w=x(3);
p=x(4);    q=x(5);     r=x(6);     om=x(4:6);
PN=x(7);   PE=x(8);    PD=x(9);    h=-PD;
e0=x(10);  e1=x(11);   e2=x(12);   e3=x(13); ev=x(10:13);
ev=ev/norm(ev); %quaternions

F=inp(27:32); % kinematic control input
V=F(1:3);     % kinematic control input (reference linear airspeed)
om=F(4:6);    % kinematic control input (reference angular airspeed)
p_c=om(1); q_c=om(2); r_c=om(3); % angular airspeed components

ws=U(1);   wh=U(2); % wind speed and heading
nu=U(3);   nv=U(4);    nw=U(5);    nn=U(3:5); % turbulence parameters

% ws,wh is incoming wind; 
wn=-ws*cos(wh); we=-ws*sin(wh); wd=0; % wind components in NED
wned=[wn we wd]';
Vw=wned;

global loc
hAGL=h-loc.hG;         % AGL altitude

%%%-------------------- Quaternions and S Transformation Matrix
[S,lam]=quaternions(e0,e1,e2,e3);

%%%-------------------- gust states and inputs (only computed if necessary)
isgust=abs(nu)>0;if isgust                     % is gust to be considered ?
   %060630 for low altitude gust is relative to wind direction  ??
   uvw_w=S*wned;ua=u-uvw_w(1);  %wind/air velocity
   % conversion from mean wind direction to body frame
   Swb=S*[cos(wh) sin(wh) 0;-sin(wh) cos(wh) 0;0 0 1];% 060630 < MIL-F-8785C 
   nn=add4noise(nn); %060925 add pseudo random 4th noise input
   [dgust_state,x_g,dx_gdt]=ddryden(x,nn,ua,hAGL);
else
   x_g=zeros(6,1);
   dgust_state=zeros(1,8);
   dx_gdt=zeros(6,1);
end
 
uvw=[u v w]';
uvw_w=S*Vw+x_g(1:3);
uvw_a=uvw-uvw_w;
ua=uvw_a(1);va=uvw_a(2);wa=uvw_a(3); % real airspeeds
%[ua(1) u uvw_w(1)-x_g(1) x_g(1) ]

% alfa and beta aerodynamic angles
Vt = sqrt( ua*ua+va*va+wa*wa ); 
if abs(ua)>.0001; aa = atan2(wa,abs(ua)); else aa=0; end;  
caa = cos(aa); saa = sin(aa);
if Vt>.0001; bb = atan2(va+1e-6,caa*ua+saa*wa); else bb=0; end;  

p=p_c;q=q_c;r=r_c-0.2*bb;
%quaternions derivative
OM4L=[0 p q r;-p 0 -r q;-q r 0 -p;-r -q p 0];
dev = -1/2*OM4L*ev + lam*ev;

%attitude angles (only for output)
tt=-asin(S(1,3)); 
ff=atan2(S(2,3),S(3,3));
yy=atan2(S(1,2),S(1,1));

V=0.05*V+0.95*uvw_a; % forgetting rate for lateral airspeed

%position derivative
dNED=S'*V+Vw; % actual ground velocities

uvw=S*dNED;   % actual local velocities
u=uvw(1); v=uvw(2); w=uvw(3);
uvwpqr=[uvw' p q r]';

ddyn_state=[uvwpqr' dNED' dev'];
other_outputs=[uvwpqr' dNED' ev' Vt [aa bb ff tt yy ua va wa]]; 

fout=[ddyn_state dgust_state other_outputs];

return

%=================
function [S,lam]=quaternions(e0,e1,e2,e3);

%outputs:
% S   transformation matrix
% lam corrective term to keep quaternion vector normalized

%quaternions and transformation matrix S
e00=e0*e0;e01=e0*e1;e02=e0*e2;e03=e0*e3;
e11=e1*e1;e12=e1*e2;e13=e1*e3;
e22=e2*e2;e23=e2*e3;
e33=e3*e3;
s11=e00+e11-e22-e33; s12=2*(e12+e03); s13=2*(e13-e02);
s21=2*(e12-e03); s22=e00-e11+e22-e33; s23=2*(e23+e01);
s31=2*(e02+e13); s32=2*(e23-e01); s33=e00-e11-e22+e33;
S=[s11 s12 s13;s21 s22 s23;s31 s32 s33];

% S is transf. matrix from fix (NED) to local frame (ABC)
lam=1-(e00+e11+e22+e33);
return

%=================
function [dgust_state,x_g,dx_gdt]=ddryden(x,nn,ua,h) 

% dryden continuous turbulence model adapted for airship case

global bm 
L=bm.L;d=bm.d;b=bm.b;

p=x(4);    q=x(5);     r=x(6);
nu=nn(1);  nv=nn(2);   nw=nn(3);

isgust=abs(nu)>0;if isgust                     % is gust to be considered ?
   xg=x(14:21); 
   Vt=max(ua,1); h=max(h,.1); % clip to avoid singularities

   xu=xg(1);   xv=xg(2:3);    xw=xg(4:5);
   xp=xg(6);    xq=xg(7);   xr=xg(8);
   
   s3=sqrt(3); 
   h0=533;  %Dryden model altitude transition (m)
   if h>h0
      lu=h0; lv=h0; lw=h0;                     %above h0=533m
   else
      lu=sqrt(h0*h); lv=lu; lw=h;
   end;
   su=sqrt(lu/h0); sv=sqrt(lv/h0); sw=sqrt(lw/h0);
   tu=lu/Vt;       tv=lv/Vt;       tw=lw/Vt;
%    tp=4*b/pi/Vt;   tq=tp*L/b;      tr=3*tq/4;  %000208:tq/tr:b->l
   tp=4*d/2/pi/Vt;   tq=tp*L/b;      tr=3*tq/4;  %l40213 para tp b->d/2
   
   dxu=-1/tu*xu+1/tu*nu*su*sqrt(2/pi*tu);
   dxv=[-2/tv -1/tv^2;1 0]*xv+[1/tv^2;0]*nv*sv*sqrt(1/pi*tv);
   dxw=[-2/tw -1/tw^2;1 0]*xw+[1/tw^2;0]*nw*sw*sqrt(1/pi*tw);
   
   yv=[s3*tv 1]*xv;
   yw=[s3*tw 1]*xw;
   
   dxp=-1/tp*xp+1/tp*nw*sw*sqrt(4/5*tp^2/tw);  %nw, not np ??? 000208
   dxq=-1/tq*xq+1/tq*(-yw);
   dxr=-1/tr*xr+1/tr*yv;
   
   ug=xu;   vg=yv;   wg=yw;
   pg=xp;
   qg=(-1/tq*xq+1/tq*(-yw))/Vt;  %sign change on 000316 /PF/JRA
   rg=(-1/tr*xr+1/tr*yv)/Vt;
   
   dgust_state=[dxu dxv' dxw' dxp dxq dxr];

%%%%%%%%%%%%%%% gust states and inputs done

   Vt1=max([1 Vt]); 
   Vg=[ug vg wg]';omg=[pg qg rg]';
   Cv=[s3*tv 1];
   Cw=[s3*tw 1];
   dVgdt=[dxu; Cv*dxv;Cw*dxw];
   domgdt =[dxp;-1/tq/Vt1*dxq-1/tq/Vt1*Cw*dxw;-1/tr/Vt1*dxr+1/tr/Vt1*Cv*dxv];
   dx_gdt=[dVgdt;domgdt];
else
   ug=0;vg=0;wg=0;pg=0;qg=0;rg=0;
   dgust_state=zeros(1,8);
   dx_gdt=zeros(6,1);
end

x_g=[ug vg wg pg qg wg]';
return

function nn=add4noise(nn)
nn4=pi*sum(nn'.*10.^(1:3)); nn4=nn4-floor(nn4); 
rms=sqrt(max(max(nn*nn')));
nn(4)=2*(nn4-.5)*rms/.83;
return