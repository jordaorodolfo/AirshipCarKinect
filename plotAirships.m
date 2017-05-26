function out = plotAirships(figureID, timePeriod, trajectory, wind, simTime)
% PLOTAIRSHIPS
% Plots the airships trajectory, using draw_droni function
% Inputs  - 
% Outputs - 
color = ['r';'b';'g'];
figure(figureID)
hold on

if timePeriod(end) > size(trajectory,1)
    out = 'Graph plot: Time period bigger than given trajectory';
    return;
end

for droniID = 1:size(trajectory,2)/3
    plotInfo(droniID) = plot(...
        trajectory(1:size(trajectory,1),(droniID-1)*3+2),...
        trajectory(1:size(trajectory,1),(droniID-1)*3+1),...
        color(droniID));
    
    for timeMoment = timePeriod
        [~,I] = min((simTime-timeMoment).^2);
        draw_droni(trajectory(I,(droniID-1)*3+1),...
            trajectory(I,(droniID-1)*3+2),...
            trajectory(I,(droniID-1)*3+3), color(droniID),0.02)
    end
end

axis equal

if wind.ws ~= 0
    axisSet = [get(gca,'xlim'),get(gca,'ylim')];
    
    [xq,yq] = meshgrid(axisSet(1)+10:20:axisSet(2),axisSet(3)+20:20:axisSet(4));

    delta_yq=-cos(wind.wh);
    delta_xq=-sin(wind.wh); 
    
    quiver(xq, yq, delta_xq*ones(size(xq)), delta_yq*ones(size(yq)),...
        'blue', 'LineWidth', 0.2, 'LineStyle',':','DisplayName','Wind',...
        'AutoScaleFactor',.9,'UData',1:10,'VData',1:10)
    text(axisSet(1)+5,axisSet(4)-5,'Wind', 'color', 'blue','FontSize', 10)
end

hold off

out = plotInfo;
end
