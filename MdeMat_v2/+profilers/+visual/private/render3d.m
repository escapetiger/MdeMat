function varAxes = render3d(axes, subdata, varName, offset)
varAxes = [];

dataNames = fieldnames(subdata);
for i = 1:length(dataNames)
    dataName = dataNames{i};
    g = subdata.(dataName);
    ax = axes(offset+i-1);
    varAxes = [varAxes, ax];

    cla(ax);
    hold(ax, 'on');

    z1 = (min(g.x1) + max(g.x1)) / 2;
    z2 = (min(g.x2) + max(g.x2)) / 2;
    z3 = (min(g.x3) + max(g.x3)) / 2;
    u = permute(g.u, [2, 1, 3]);
    h = slice(ax, g.x1, g.x2, g.x3, u, z1, z2, z3);
    set(h, 'EdgeColor', 'none');

    xlabel(ax, 'x');
    ylabel(ax, 'y');
    zlabel(ax, 'z');
    colormap(ax, 'turbo');
    colorbar(ax);
    axis(ax, 'vis3d', 'equal', 'tight');
    setColorLimits(ax, u);
    title(ax, sprintf('%s: %s', strrep(dataName, '_', '-'), varName));
    view(ax, 3);

    hold(ax, 'off');
end
end
