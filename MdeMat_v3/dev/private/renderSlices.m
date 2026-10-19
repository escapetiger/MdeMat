function renderSlices(axes, subdata, styleset, varName, offset, nDims)
dataNames = fieldnames(subdata);
if length(dataNames) < 1
    return;
end

ax = axes(offset+length(dataNames));
cla(ax);
hold(ax, 'on');

[xmin, xmax, ymin, ymax] = deal(inf, -inf, inf, -inf);
legends = cell(1, length(subdata));

for i = 1:length(dataNames)
    dataName = dataNames{i};
    g = subdata.(dataName);

    r = linspace(min(g.x1), max(g.x1), 64);
    if nDims == 2
        R = sqrt(2);
        s = interp2(g.x1, g.x2, g.u.', r/R, r/R);
    elseif nDims == 3
        R = sqrt(3);
        u = permute(g.u, [2, 1, 3]);
        s = interp3(g.x1, g.x2, g.x3, u, r/R, r/R, r/R);
    end

    if isfield(styleset, dataName)
        style = styleset.(dataName);
    else
        style = {};
    end

    plot(ax, r, s, style{:}, 'DisplayName', dataName);

    xmin = min(xmin, min(r(:)));
    xmax = max(xmax, max(r(:)));
    if ~isempty(s)
        ymin = min(ymin, min(s(:)));
        ymax = max(ymax, max(s(:)));
    end

    legends{i} = strrep(dataName, '_', '-');
end

if ymin == ymax
    ymin = ymin - 1e-8;
    ymax = ymax + 1e-8;
end
xlim(ax, [xmin, xmax]);
ylim(ax, [ymin, ymax]);
xlabel(ax, 'x');
ylabel(ax, varName);
legend(ax, legends, 'Location', 'best');

title(ax, 'diagonal slice');

hold(ax, 'off');
end