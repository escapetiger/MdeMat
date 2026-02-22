function ax = render1d(axes, subdata, styleset, varName, offset)
dataNames = fieldnames(subdata);
ax = axes(offset);
cla(ax);
hold(ax, 'on');

[xmin, xmax, umin, umax] = deal(inf, -inf, inf, -inf);
legends = cell(1, length(dataNames));

for i = 1:length(dataNames)
    dataName = dataNames{i};
    g = subdata.(dataName);

    xmin = min(xmin, min(g.x1));
    xmax = max(xmax, max(g.x1));
    umin = min(umin, min(g.u));
    umax = max(umax, max(g.u));

    if isfield(styleset, dataName)
        style = styleset.(dataName);
    else
        style = {};
    end

    plot(ax, g.x1(:).', g.u(:).', style{:});
    legends{i} = strrep(dataName, '_', '-');
end

xlim(ax, [xmin, xmax]);
if umax > umin
    ylim(ax, [umin, umax]);
end
xlabel(ax, 'x');
ylabel(ax, varName);
legend(ax, legends, 'Location', 'best');

hold(ax, 'off');

end

