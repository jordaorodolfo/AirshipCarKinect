function out = plotMovie(figureID, timePeriod, trajectory, wind, simTime, axisSet)
% PLOTAIRSHIPS
% Plots the airships trajectory, using draw_droni function
% Inputs  -
% Outputs -
color = ['r';'b';'g'];
figure(figureID)

if timePeriod(end) > size(trajectory,1)
    out = 'Graph plot: Time period bigger than given trajectory';
    return;
end

n = 1;
for timeMoment = timePeriod
    [~,I] = min((simTime-timeMoment).^2);
    
    clf
    hold on
    
    plot(trajectory(1:I,2), trajectory(1:I,1), 'r');
    plot(trajectory(1:I,5), trajectory(1:I,4), 'b');
    plot(trajectory(1:I,8), trajectory(1:I,7), 'g');
    plot(trajectory(1:I,11), trajectory(1:I,10), 'black-.');

    axis(axisSet)
    
    draw_droni(trajectory(I,1),...
        trajectory(I,2),...
        trajectory(I,3), color(1),0.02)
    
    draw_droni(trajectory(I,4),...
        trajectory(I,5),...
        trajectory(I,6), color(2),0.02)
    
    draw_droni(trajectory(I,7),...
        trajectory(I,8),...
        trajectory(I,9), color(3),0.02)
    
    [xq,yq] = meshgrid((-200:20:500)-timeMoment*sind(30),...
        (-200:20:200)-timeMoment*cosd(30));
    
    delta_yq=-cos(wind.wh);
    delta_xq=-sin(wind.wh);
    
    quiver(xq, yq, delta_xq*ones(size(xq)), delta_yq*ones(size(yq)),...
        'blue', 'LineWidth', 0.2, 'LineStyle',':','DisplayName','Wind')
    text(-150,40,'Wind', 'color', 'blue','FontSize', 10)

    hold off
    
    movieFrame(n) = getframe;
    n = n+1;
end

myVideo = VideoWriter('movingTarget.avi');
myVideo.FrameRate = 7.5;  % Default 30
myVideo.Quality = 90;    % Default 75
open(myVideo);
writeVideo(myVideo, movieFrame);
close(myVideo);

out = movieFrame;
end