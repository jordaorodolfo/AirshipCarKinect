function [nn,dd,ee,omdt,imis,Pr,llmis,smis]=misptrk(N,E,GS,ll,Pr0,ds,mis,imis,verbose)

%--function [nn,dd,ee,om,imis,Pr,llmis]=misptrk(N,E,GS,ll,Pr0,ds,mis,imis,verbose)
%
%-- mission horiz. path tracking errors (nn,dd,ee) and reference yaw rate (om)
%-- + updated mission index (imis) and current tracked point (Pr)
%
%inputs:
%-- (N,E) measured position
%-- (GS,ll) measured ground speed and heading
%-- Pr0=(Nr0,Er0) previous tracked point
%-- (ds) curvilinear length to next tracked point Pr=(Nr,Er)
%-- mis=[N E u h du dh tt t] mission
%-- (imis) index of current mission segment
%
%outputs:
%-- (nn,dd) longitudinal and lateral position tracking error
%-- (ee) angular tracking error
%-- (om) desired yaw rate
%-- (imis) updated mission index
%-- Pr=(Nr,Er) updated tracked point
%
%040303:JRA update from mistrk for path tracking w/longitudinal error nn

P=[N E]';Pr0=Pr0(:);
ext=[0 -1;1 0];  %matrix 90 deg rotation for cross product 040303
A=mis(imis,1:2)'; B=mis(imis+1,1:2)';
ab=B-A; dab=sqrt(ab'*ab);               %straight segment length
pb=B-P; dpb=sqrt(pb'*pb);               %dist. to 2nd point (b)
if dab==0; A=P;ab=pb;dab=dpb;end        %040315-----to cope w/ hovering 
ttab=mis(imis,7); 
if (ttab==0)                            % straight segment case
  Pr=Pr0+ds*ab/dab;
  bPr=Pr-B;dbPr=sqrt(bPr'*bPr);         %dist. from tracked point to (b)
  cb=bPr'*(-ab);                        %dist. along segment (cos)
else
  M=(A+B)./2;
  C=[M(1);M(2)]+ext*[ab(1);ab(2)]./tan(ttab/2)./2;
  sgn=(2*(ttab>=0)-1);
  r=dab/2/sin(ttab/2); %signed radius
  ttr=ds/r;Rttr=[cos(ttr) sin(ttr);-sin(ttr) cos(ttr)]';
  Pr=C+Rttr*(Pr0-C);
  bPr=Pr-B;dbPr=sqrt(bPr'*bPr);          %dist. from tracked point to (b)
  cb=bPr'*(Pr0-B);                      %dist. along segment (cos)
%   global XX
  CA=A-C;CPr=Pr-C;cb=abs(ttab)-abs(atan2(-CA'*ext*CPr,CA'*CPr));%keyboard
%   XX=[XX;imis atan2(-CA'*ext*CPr,CA'*CPr) cb];
end
%   if imis==2,keyboard,end

if (cb<0);               %then segment is over, next one
% cbmin=dab/50;if (cb<cbmin);  %120403: allow interval for seg.end: implies some step
  [np,mp]=size(mis);
  if imis == np-1
    imis=1;                            %for end of mission, last point is first
  else                                 %else switch to next segment
    imis=imis+1;
    A=mis(imis,1:2)';   B=mis(imis+1,1:2)';
    ab=B-A; dab=sqrt(ab'*ab);
    pb=B-P; dpb=sqrt(pb'*pb);
if dab==0; A=P;ab=pb;dab=dpb;end        %040315-----to cope w/ hovering 
    ttab=mis(imis,7); 
    if (nargin>8)&(verbose); 
       st=sprintf('+++> misptrk: mis=[%6.0f %6.0f %6.0f %6.0f %6.0f] dpb=%6.1fm',...
          mis(imis,[1:4 7]),dpb);disp(st),
    end
  end;
  ds1=dbPr;
  if (ttab==0)                            % straight segment case
    Pr=A+ds1*ab/dab;
  else
    M=(A+B)./2;
    C=[M(1);M(2)]+ext*[ab(1);ab(2)]./tan(ttab/2)./2;
    sgn=(2*(ttab>=0)-1);
    r=dab/2/sin(ttab/2); %signed radius
    ttr=ds1/r;Rttr=[cos(ttr) sin(ttr);-sin(ttr) cos(ttr)]';
    Pr=C+Rttr*(A-C);
  end
end;

%%%%%--------------------------------------tracking error
p0=P-GS*[cos(ll);sin(ll)];                 % position 1s before
p0p=P-p0;dp=sqrt(p0p'*p0p);  
if (ttab==0)                               % straight segment case
  %dd0=ab'*ext*pb/dab;                     % dist. to segment
  t=ab/dab;                                % unitary tangent vector
  omdt=0;
  smis=mislen(mis(1:imis,:))+sqrt((Pr-A)'*(Pr-A));
  %disp(1);keyboard
else                                       % circular segment case
  %cp=P-C;dd0=-sgn*(sqrt(cp'*cp)-r*sgn);   % distance to circle segment
  t=ext*(Pr-C)*sgn/abs(r);                 % unitary tangent vector
  omdt=ds/r; %keyboard
  CA=A-C;CPr=Pr-C;ttAP=abs(atan2(-CA'*ext*CPr,CA'*CPr));%keyboard
  smis=mislen(mis(1:imis,:))+ttAP*abs(r);
end

llmis=atan2(t(2),t(1));

% position and angular tracking errors
nn=(P-Pr)'*t;
dd=(P-Pr)'*ext*t;
see=p0p'*ext*t;  cee=p0p'*t;
ee=atan2(see,cee);

%keyboard 