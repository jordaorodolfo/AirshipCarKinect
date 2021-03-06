function [Xge]=ftrack2(x_y_yy, t)
    %*function [xe, ye, yye]=ftrack_adapt(x_y_yy,t)
    %-- 2-D mission path tracking errors (xe,ye,yye))
    %-- updated mission index (imis) and current tracked point (Pr)
    %
    %inputs:
    %-- x_y_yy: Vector with 2-D position and yaw angle of the vehicle
    %-- t: Simulation time (use simulink clock block)
    %
    %outputs:
    %-- (nn,dd) longitudinal and lateral position tracking error
    %-- (ee) angular tracking error
    %
    %global variables
    %-- Pr0: Previous tracked position
    %-- t0: Previous sampled time
    %-- uref: reference speed (need to be defined in main program)
    %-- mis: Mission matrix (np x 3) with np points (xy position) and arc travelled (0 if straight line)
    %-- imis: actual mission point index
    %
    global Pr0 t0 uref mis imis
    
    if t==0                                %First iteration
        t0=0;
        Pr0=mis(1,1:2);
        imis=1;
    end

    % Initialize the function
    P=zeros(2,1);
    P(1)=x_y_yy(1);P(2)=x_y_yy(2);          %Vehicle 2-D position
    yy=x_y_yy(3);                           %Vehicle yaw
    ds=uref*(t-t0);                         %Distance travelled by the reference vehicle
    ext=[0 -1;1 0];                         %90deg rotation matrix for cross product
    A=mis(imis,1:2)'; B=mis(imis+1,1:2)';   %Start and finish of the mission segment
    ab=B-A; dab=sqrt(ab'*ab);               %straight segment length
    pb=B-P; dpb=sqrt(pb'*pb);               %distance to finish point (B)
    ttab=mis(imis,3);                       %angular variation of the track (0 if straight line)

    %Verifying if the reference is still in the previous track segment
    if (ttab==0)                            %straight segment case
      Pr=Pr0+ds*ab/dab;                     %Tracked point (previous track segment)
      bPr=Pr-B;dbPr=sqrt(bPr'*bPr);         %distance from tracked point to (B)
      cb=bPr'*(-ab);                        %distance along segment (cos)
    else                                    %circular arc case
      M=(A+B)./2;
      C=[M(1);M(2)]+ext*[ab(1);ab(2)]./tan(ttab/2)./2;  %Center of the circle
      sgn=(2*(ttab>=0)-1);
      r=dab/2/sin(ttab/2);                  %signed radius
      ttr=ds/r;                             %arc travelled
      Rttr=[cos(ttr) sin(ttr);-sin(ttr) cos(ttr)]'; %Rotation matrix of the travelled arc
      Pr=C+Rttr*(Pr0-C);                    %Tracked point (previous track segment)
      bPr=Pr-B;dbPr=sqrt(bPr'*bPr);         %distance from tracked point to (B)
      CA=A-C;CPr=Pr-C;
      cb=abs(ttab)-abs(atan2(-CA'*ext*CPr,CA'*CPr));%difference between arcs
    end

    if (cb<0);                              %then segment is over, next one
      [np,mp]=size(mis);
      if imis == np-1
        imis=1;                             %for end of mission, last point is first
      else                                  %else switch to next segment
        imis=imis+1;
        A=mis(imis,1:2)';   B=mis(imis+1,1:2)';
        ab=B-A; dab=sqrt(ab'*ab);
        pb=B-P; dpb=sqrt(pb'*pb);
        ttab=mis(imis,3); 
      end;
      ds1=dbPr;                             %part of the distance travelled in the actual track segment
      if (ttab==0)                          %straight segment case
        Pr=A+ds1*ab/dab;                    %Tracked point (new track segment)
      else                                  %circular arc case
        M=(A+B)./2;
        C=[M(1);M(2)]+ext*[ab(1);ab(2)]./tan(ttab/2)./2; %Center of the circle
        sgn=(2*(ttab>=0)-1);
        r=dab/2/sin(ttab/2);                %signed radius
        ttr=ds1/r;                          %arc travelled
        Rttr=[cos(ttr) sin(ttr);-sin(ttr) cos(ttr)]'; %Rotation matrix of the travelled arc
        Pr=C+Rttr*(A-C);                    %Tracked point (new track segment)
      end
    end;

    %Calculating the tracking error
    p0=P-uref*[cos(yy);sin(yy)];            %Position 1seg before
    p0p=P-p0;dp=sqrt(p0p'*p0p);             %Vector with the corresponding yaw
    %Obtaining the tangent vector of the tracked point
    if (ttab==0)                            %straight segment case
      t=ab/dab;                             %unitary tangent vector
    else                                    %circular segment case
      t=ext*(Pr-C)*sgn/abs(r);              %unitary tangent vector
    end
    % position and angular tracking errors
    xe=(P-Pr)'*t;                           %Longitudinal projection (xe)
    ye=(P-Pr)'*ext*t;                       %Lateral projection (ye)                      
    yproj_yy=p0p'*ext*t;  xproj_yy=p0p'*t;  %Vehicle vector projections
    yye=atan2(yproj_yy,xproj_yy);           %Angle between track tangent and vehicle
    % Update global variables
    Pr0=Pr;
    t0=t;
 
 % Controle PID n�o linear   
cs=1/6 ; % inverso do raio
Kp=0.5;Kd=1.2;
L=1.3;
cy1=1-cs*ye;
cce=cos(yye);
tte=tan(yye);
delta=atan(L*((cce^3/cy1)*(-Kd*cy1*tte-Kp*ye+cs*cy1*tte*tte)+cs*cce/cy1));

% a2=ye;
% a3=cy1*tte;
% z=Lambda*a2+a3;
% %delta=atan(L*((cce^3/cy1^2)*(-K*z-0*Lambda*a3-rho*tanh(0.5*z)+cs*cy1*tte*tte)+cs*cce/cy1));

Xge=[delta,xe,ye,yye]';
end