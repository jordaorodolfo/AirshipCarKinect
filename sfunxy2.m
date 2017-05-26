function  [sys, x0]  = sfunxy2(~,x,u,flag,ax,autoscale)
%SFUNXY an S-function which acts an X-Y scope using MATLAB's plotting functions.
%
%	This M-file is designed to be used in a SIMULINK S-function block.
%	It stores the last input point using discrete states 
%	and then plots this against the most recent input. 
%
%	Set this M-file up in an S-function block with a MUX with two
%	inputs feeding into this block. Set the function parameter
%	up as a four element vector which defines the axis of the graph.
%
%	See also SFUNXYS, LORENZS.

%	Copyright (c) 1991-93 by the MathWorks, Inc.
%	Andrew Grace 05-30-91.
%   Pedro Gatti  03-06-17

if abs(flag) == 2   	% A real time hit 
	if (x(1) ~= Inf) 
		% Use none as Erasemode for fast plotting
		plot([u(2),x(2)],[u(1),x(1)],'k-');
        plot([u(4),x(4)],[u(3),x(3)],'r-');
        plot([u(6),x(6)],[u(5),x(5)],'b-');
        plot([u(8),x(8)],[u(7),x(7)],'g-');
        if autoscale == 1
            axis auto		% Sets autoscale on
        end
	end
	sys = u;
elseif flag  == 0 	% Initialization
        sys = [0;8;0;8;0;0]; 	% System sizes - 2 discrete states, 2 inputs
        x0 = [Inf;Inf;Inf;Inf;Inf;Inf;Inf;Inf]; 	% Flag to indicate first point
	% Initialize graph;
	hold off
	figure(1);clf
    set(gcf,'Position',[500+350 375 476 441]); %[850 376 476 441]
	if nargin < 5  
		disp('Axis not defined; using current axis');
		ax = axis;
	elseif length(ax)~=4
		disp('Axis not defined correctly; it must be a 1x4 vector')
		ax = axis;    
	end
	axis(ax)		% Set up the axis
    axis square
	hold on
	plot(ax(1),ax(3))  	% Force a box to be drawn
	xlabel('East [m]');	ylabel('North [m]');
    %plot(xy(2,:),xy(1,:));
    %plot(xy(1,:),xy(2,:),'g-',xy(3,:),xy(4,:),'r-',xy(5,:),xy(6,:),'b-')
    set(gcf,'Name','East-North')
else 
        sys = [];
end
