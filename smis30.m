function [sys,x0]=smis30(t,x,U,flag,param,mis,dt0)

%function [sys,x0]=smis30(t,x,U,flag,param,mis,dt0)
%
% mission error, tracking and references for sbm26
%
%990630: circular segments pre-computed
%990615: GS output added
%991116: dte computed and output
%000214: teps0 instead of eps0 (eps0=GS*teps0)
%
% state is x=[t0;p0x;p0y;imis;dd;ee;cx;cy;sgn;r;om;GS;ucmd;hcmd;dte] (14+1)
% input is U=[u;h;N;E] (4)
% output   y=[imis ucmd hcmd dd ee om GS t0 +dte] (8+1)
%      param=[uref Ni Ei hi teps0 ui]
%        mis=[N E u h du dh tt t]

% next discrete state %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if abs(flag)==2
    t0=x(1);
    if ((t-t0)>=dt0)		%then ok new step comp. (simulink step time control)
        %------------------------------------------------
        dt=t-t0;
        p0=[x(2);x(3)];imis=x(4);
        c=x(7:8);
        psi_traj=x(7);
        sgn=x(9);r=x(10);GS=x(12);
        u0=x(13);h0=x(14);
        u=U(1);h=U(2);p=U(3:4);
        teps0=param(5);
        
        %  p0
        %  p
        
        %eps0=max(GS,.1)*teps0;
        eps0=10;
        
        
        %%%%%--------------------------------------mission control
        ext=[0 -1;1 0];
        a=mis(imis,1:2)'; b=mis(imis+1,1:2)';
        ab=b-a; dab=sqrt(ab'*ab);               %straight segment length
        pb=b-p; sb=sqrt(pb'*pb);                %dist. to 2nd point (b)
        cb=pb'*ab/dab/sb;                      %dist. along segment (cos)
        %a
        %b
        
        if (sb<eps0) | (cb<0);
            [np,mp]=size(mis);
            if imis == np-1
                imis=1;                           %for end of mission, last point is first
            else                                 %else switch to next segment
                imis=imis+1;
                
                a=mis(imis,1:2)';   b=mis(imis+1,1:2)';
                ab=b-a; dab=sqrt(ab'*ab);
                pb=b-p;
                
                straj=ab'*ext*[1;0];
                ctraj=ab'*[1;0];
                psi_traj=atan2(straj,ctraj);
                psi_traj=psi_traj*180/pi;
            end
        end
        
        %%% circular segment pre-computing
        ttab=mis(imis,7);
        if (abs(ttab)>0)
            m=(a+b)./2;
            c=[m(1);m(2)]+ext*[ab(1);ab(2)]./tan(ttab/2)./2;
            sgn=2*(ttab>=0)-1;
            r=dab/2/sin(ttab/2);
        else
            c=[0;0];sgn=1;r=1;
        end;
        %disp(mis(imis,1:4));
        
        %%%%%--------------------------------------lateral tracking and error
        p0p=p-p0;   dp=sqrt(p0p'*p0p);  GS=dp/dt;
        dd=pb'*ext*ab/dab;                      % dist. to segment
        %dd
        
        if (dp<dt0*.5)                          % should not occur
            ee=0;  om=0; GS=x(12);
        else                                    % straight segment case
            see=ab'*ext*p0p;
            cee=ab'*p0p;
            ee=atan2(see,cee);
            ttab=mis(imis,7);om=0;
            % circular segment case
            if (abs(ttab)>0)
                % p
                % c
                % r
                cp=p-c;
                dd=(cp'*cp-r^2)/2/r;               % linearized distance error
                
                cee=p0p'*cp*sgn;  see=p0p'*ext*cp*sgn;
                ee=atan2(cee,see);
                om=GS/r;
            else Tom=.3;om=x(11)*(1-dt/Tom);%om=0;     % low-pass to quit circle
            end;
        end;
        %%%%%--------------------------------------longitudinal command
        
        umin=4;    umax=13;         % air speed interval for lon. stability (bmlon20h)
        dumax=.5;  dhmax=1;         % maximum lon. reference rates
        
        V=mis(imis,3);  H=mis(imis,4);
        if V>umax; V=umax; elseif V<umin;  V=umin;  end;
        DV=mis(imis,5); DH=mis(imis,6);
        DV=abs(DV)*(2*(V>u0)-1);DH=abs(DH)*(2*(H>h0)-1);  %guarantee correct sign rates
        if DV>dumax;DV=dumax;elseif DV<-dumax;DV=-dumax;end;
        if DH>dhmax;DH=dhmax;elseif DH<-dhmax;DH=-dhmax;end;
        
        if u<.5; DH=0;  end;
        
        if abs(V-u0)<abs(DV*dt);ucmd=V;
        else ucmd=u0+dt*DV;
        end;
        if abs(H-h0)<abs(DH*dt);hcmd=H;
        else hcmd=h0+dt*DH;
        end;
        
        dte=sb;[nmis,mmis]=size(mis);
        for ii=imis+1:nmis-1;
            a=mis(ii,1:2)';   b=mis(ii+1,1:2)';
            ab=b-a; dab=sqrt(ab'*ab);
            dte=dte+dab;
        end;
        dte0=x(15);dte=0.5*(dte+dte0);%if abs(dte0-dte)>3*GS*dt; dte=dte0; end;
        if dte<eps0*1.1; disp('mission over.');end;
        
        %%%%%--------------------------------------
        % state is x=[t0;p0x;p0y;imis;dd;ee;cx;cy;sgn;r;om;GS;ucmd;hcmd] (14)
        
        % sys=[t;p(1);p(2);imis;dd;ee;c(1);c(2);sgn;r;om;GS;ucmd;hcmd;dte];
        sys=[t;p(1);p(2);imis;dd;ee;psi_traj;c(2);sgn;r;om;GS;ucmd;hcmd;dte];
        
    else
        sys=x;
    end;
    
    % output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif flag==3
    
    t0=x(1);    imis=x(4);                         % control
    ucmd=x(13);
    hcmd=x(14);                     % lon ref.
    dd=x(5);    ee=x(6);    om=x(11);    GS=x(12); % lat.outputs
    dte=x(15);
    psi_traj=x(7);
    
    %sys=[imis ucmd hcmd dd ee om GS t0 dte];
    sys=[psi_traj ucmd hcmd dd ee om GS t0 dte];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif flag==4
    t0=floor(t/dt0+1e-9)*dt0;
    sys=t0+dt0;
    
    %init sizes and values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif flag==0
    %  sizes: #cont st; #disc.st.; #outputs; #inputs; 0;1
    
    sys=[0 14+1 8+1 4,0,1];
    
    %initial state value
    % state is x=[t0;p0x;p0y;imis;dd;ee;cx;cy;sgn;r;om;GS;ucmd;hcmd;dte] (14+1)
    imis=1;
    dte=0;[nmis,mmis]=size(mis);
    for ii=imis:nmis-1;
        a=mis(ii,1:2)';   b=mis(ii+1,1:2)';
        ab=b-a; dab=sqrt(ab'*ab);  %as if straight lines only
        dte=dte+dab;
    end;
    
    Ni=param(2);Ei=param(3);hi=param(4);ui=param(6);
    p0x=Ni;p0y=Ei;ucmd=ui;hcmd=hi;
    
    
    x0=[0 p0x p0y imis 0 0 0 1 1 0 0 ucmd ucmd hcmd dte];
    
end;  %of all
