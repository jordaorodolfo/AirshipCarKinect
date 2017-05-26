function [] = draw_droni(x, y, yy, color, scale)
%DRAW_DRONI Draw an airship in the given position (x,y) with the given yaw
%angle (yy) and the given scale.
%
%   x - North position
%   y - East position
%   yy - yaw angle (rad)
%   scale - airship scale (smaler scale => smaller airship), usually
%   between 0 and 1.
%%%
delta_x = 300*scale;
delta_y = 600*scale;

yy=atan2(sin(yy+pi),cos(yy+pi));

%% Envelope
vals = [218.57143,588.07649 4.28571,-200 84.28571,-201.42857 84.28571,-201.42857;
386.60151,588.07649 -4.28571,-200 -84.28571,-201.42857 -84.28571,-201.42857;
386.60151,588.07649 1.01016,180.81731 -60.10407,296.98485 -60.10407,296.98485;
218.57303,588.07649 -1.01016,180.81731 60.10407,296.98485 60.10407,296.98485;
278.60007,884.68951 14.40731,23.68808 25.01391,21.2132 25.01391,21.2132;
326.57074,884.68951 -14.40731,23.68808 -25.01391,21.2132 -25.01391,21.2132];

vals = vals*scale;
vals(:,1) = vals(:,1)-delta_x;
vals(:,2) = vals(:,2)-delta_y;

for idx = 1:size(vals,1)
    p = [[vals(idx,1),vals(idx,3:2:end)+vals(idx,1)]',[vals(idx,2),vals(idx,4:2:end)+vals(idx,2)]'];
    n = size(p,1);
    n1=n-1;
    for  i=0:1:n1
        sigma(i+1)=factorial(n1)/(factorial(i)*factorial(n1-i));  % for calculating (x!/(y!(x-y)!)) values 
    end
    l=[];
    UB=[];
    for u=0:0.002:1
        for d=1:n
            UB(d)=sigma(d)*((1-u)^(n-d))*(u^(d-1));
        end
        l=cat(1,l,UB);                                      %catenation 
    end
    P=l*p;
    hl = line(P(:,1)*cos(yy)+P(:,2)*sin(yy)+y,P(:,2)*cos(yy)-P(:,1)*sin(yy)+x, 'color', color, 'LineWidth', 1);
end
%%%
%% Calda lateral
vals = [238.00839,764.57617 -26.92692,88.21712 55.64897,1.98688;
367.67063,764.57617 26.92692,88.21712 -55.64897,1.98688];

vals = vals*scale;
vals(:,1) = vals(:,1)-delta_x;
vals(:,2) = vals(:,2)-delta_y;

for idx=1:size(vals,1)

P = [cumsum(vals(idx,1:2:end))',cumsum(vals(idx,2:2:end))'];
hl = plot(P(:,1)*cos(yy)+P(:,2)*sin(yy)+y,P(:,2)*cos(yy)-P(:,1)*sin(yy)+x, 'color', color, 'LineWidth', 1);

end
%%%
%% Calda vertical
vals=[300.75,855.36221 0,-99.25];

vals = vals*scale;
vals(:,1) = vals(:,1)-delta_x;
vals(:,2) = vals(:,2)-delta_y;

P = [cumsum(vals(1:2:end))',cumsum(vals(2:2:end))'];
hl = plot(P(:,1)*cos(yy)+P(:,2)*sin(yy)+y,P(:,2)*cos(yy)-P(:,1)*sin(yy)+x, 'color', color, 'LineWidth', 1);

%%%
end

