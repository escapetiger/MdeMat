function varAxes = render2d(axes, subdata, varName, offset)
dataNames = fieldnames(subdata);
varAxes = [];
for i = 1:length(dataNames)
    dataName = dataNames{i};
    g = subdata.(dataName);
    ax = axes(offset+i-1);
    varAxes = [varAxes, ax];

    cla(ax);
    hold(ax, 'on');

    imagesc(ax, g.x1, g.x2, g.u');
    axis(ax, 'xy', 'equal', 'tight');
    xlabel(ax, 'x');
    ylabel(ax, 'y');
    colormap(ax, 'jet');
    colorbar(ax);

    setColorLimits(ax, g.u);

    title(ax, sprintf('%s (%s)', varName, dataName));

    hold(ax, 'off');
end
end
