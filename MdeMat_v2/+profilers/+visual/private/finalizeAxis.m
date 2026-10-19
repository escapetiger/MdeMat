function finalizeAxis(ax, xmin, xmax, ymin, ymax, xLabel, yLabel, legends)
% FINALIZEAXIS Applies common formatting to an axis.
%
% Syntax:
%   finalize(ax, xmin, xmax, ymin, ymax, xLabel, yLabel, legends)
%
% Inputs:
%   ax - Axes handle to format
%   xmin - Minimum x-axis value
%   xmax - Maximum x-axis value
%   ymin - Minimum y-axis value
%   ymax - Maximum y-axis value
%   xLabel - Label for x-axis
%   yLabel - Label for y-axis
%   legends - Cell array of legend entries

epsilon = 1e-8;
xmin = min(-epsilon, xmin);
xmax = max(epsilon, xmax);
ymin = min(-epsilon, ymin);
ymax = max(epsilon, ymax);

xlim(ax, [xmin, xmax]);
ylim(ax, [ymin, ymax]);
xlabel(ax, xLabel);
ylabel(ax, yLabel);

if ~isempty(legends)
    legend(ax, legends, 'Location', 'best');
end

hold(ax, 'off');
end