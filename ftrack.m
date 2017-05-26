function Xge=ftrack(p)

% mission tracking error to test car guidance for fdlat/sdlatxy
%090925

global mis imis dmis circ 

Kp=0.5;Kd=1.2;
%Kp=0.1;Kd=0.8;
L=1.3;

%Lambda=0.1;
Lambda=0.4167;
%K=0.2; 
K=1.2;
rho=0.45;

xy=p(1:2);psi=p(3);

nmis=size(mis,1);
b=mis(imis+1,:)';a=mis(imis,:)';
dpb=norm(b-xy);
if dpb<dmis; 
    imis=imis+1;
    if imis==nmis;imis=1;end
end
ab=b-a;dab=norm(ab);pb=b-xy;
ye=[0 0 1]*cross([pb;0],[ab;0])/dab;
psie=psi-atan2(ab(2),ab(1));
cs=0;

   if psie>5       
        psie=psie-2*pi;     
    elseif psie<-5
        psie=psie+2*pi;  
    end

% if circ &(rem(imis,2)==0)
%     % lower circle not fully ok 090925
%     if imis==2;c=[100 50]';else c=[0 50]';end
%     cp=xy-c;dcp=norm(cp);
%     ye=50-dcp;
%     theta=atan2(cp(2),cp(1));
%     psie=psi-theta-pi/2;
% end
theta=0;
if circ &(rem(imis,2)==0)
    % lower circle not fully ok 090925
   
    raio=6;
    if imis==2;c=[10 raio]';else c=[0 raio]';end
    cp=xy-c;dcp=norm(cp);
    ye=raio-dcp;
    theta=atan2(cp(2),cp(1));
         
    psie=psi-theta-pi/2;
    
     if psie>5       
        psie=psie-2*pi;     
    elseif psie<-5
        psie=psie+2*pi;  
    end
          
        raio=6;
    cs=1/raio;
end

cy1=1-cs*ye;
cce=cos(psie);
tte=tan(psie);
delta=atan(L*((cce^3/cy1)*(-Kd*cy1*tte-Kp*ye+cs*cy1*tte*tte)+cs*cce/cy1));

a2=ye;
a3=cy1*tte;
z=Lambda*a2+a3;
%delta=atan(L*((cce^3/cy1^2)*(-K*z-0*Lambda*a3-rho*tanh(0.5*z)+cs*cy1*tte*tte)+cs*cce/cy1));

Xge=[delta ye psie]';
        