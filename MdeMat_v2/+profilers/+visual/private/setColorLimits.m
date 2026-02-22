function setColorLimits(ax, u)
% SETCOLORLIMITS Sets color axis limits based on data values.
%
% Syntax:
%   setColorLimits(ax, u)
%
% Inputs:
%   ax - Axes handle to set color limits on
%   u - Data array used to determine limits

uMin = min(u(:));
uMax = max(u(:));
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

caxis(ax, cLimits);
end
