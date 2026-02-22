function syncColorLimits(axes, subdata)
% SYNCCOLORLIMITS Sets unified color limits for multiple plots.
%
% Syntax:
%   syncColorLimits(axes, subdata)
%
% Inputs:
%   axes - Array of axes handles to synchronize
%   subdata - Map of data names to structures containing data

if isempty(fieldnames(subdata))
    return;
end

uMin = inf;
uMax = -inf;

dataNames = fieldnames(subdata);
for i = 1:length(dataNames)
    dataName = dataNames{i};
    u = subdata.(dataName).u;
    uMin = min(uMin, min(u(:)));
    uMax = max(uMax, max(u(:)));
end

absMax = max(abs(uMin), abs(uMax));

if uMin < 0 && uMax > 0
    cLimits = [-absMax, absMax];
else
    padding = 0.05 * (uMax - uMin);
    if padding == 0
        padding = 0.1 * abs(uMin);
        if padding == 0
            padding = 0.1;
        end
    end
    cLimits = [uMin - padding, uMax + padding];
end

for i = 1:length(axes)
    caxis(axes(i), cLimits);
end
end
