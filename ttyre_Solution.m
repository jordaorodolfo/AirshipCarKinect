if model==1 %Pacejka Model

aab=Cy*tan(aa)*15; %corrected side slip angle for magick formula

%combined model variables
a1=aab/aam;s=sig/sm;rho=sqrt(a1.^2+s.^2);rho=max(rho,1e-9*[1 1 1 1]);
sx=rho*sxm;aab=rho*aam;
Dx=mu*Fz;
Fx1=pacejka(sx,0,Bx,Cx,Dx,Ex,0);
Dy=mu*Fz;
Fy1=pacejka(aab,0,By,Cy,Dy,Ey,0);
Dr=Dy*aw/3; %101019:JRA Rajamani p420/ aw = half length of contact patch
Mz1=pacejka(aab,0,Br,Cr,Dr,Er,0);

Fx=Fx1.*s./rho;Fy=Fy1.*a1./rho;Mz=Mz1.*a1./rho;

else  % Dugoff Model
    
    for i=1:4
        muFz=mu*Fz(i);      
        Cs=100*muFz*Bx*Cx;
        Ca=muFz*By*Cy;  
       
        Fr=sqrt(((Cs*sig(i))^2)+(Ca*tan(aa(i)))^2);
                        
        if Fr>muFz/2
            f=muFz/(2*Fr)*(2-muFz/(2*Fr)); 
        else 
            f=1; 
        end 
        
        Fy(i)=Ca*tan(aa(i))*f/(1+sig(i));
        Fx(i)=Cs*sig(i)*f/(1+sig(i));        
    end
end
return

function Y=pacejka(X,Sx,B,C,D,E,Sy)
% Pacejka magic formula
x=X+Sx;
Y=D.*sin( C*atan( B*x - E*( B*x - atan(B*x)))) +Sy;
return