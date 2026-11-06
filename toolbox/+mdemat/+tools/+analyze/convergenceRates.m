function report = convergenceRates(scales, errors)
%CONVERGENCERATES Compute observed rates for an explicitly supplied study.
%   Scales H is positive and strictly decreasing. Errors E is a nonnegative
%   finite matrix [refinements, fields], or one vector for a single field.
%   LocalRates(i,j) = log(E(i,j)/E(i+1,j))/log(H(i)/H(i+1)).
%   A rate involving an exactly zero error is NaN and Defined is false.
%   Negative rates remain visible. This function does not run a solver or
%   infer which discretization/temporal/rank error the caller measured.
arguments
    scales {mustBeNumeric, mustBeReal, mustBeVector, mustBePositive, mustBeFinite}
    errors {mustBeNumeric, mustBeReal, mustBeNonnegative, mustBeFinite}
end
scales = double(scales(:));
if isvector(errors), errors = errors(:); end
errors = double(errors);
if numel(scales) < 2 || any(diff(scales) >= 0) || ...
        ~ismatrix(errors) || size(errors, 1) ~= numel(scales) || size(errors, 2) == 0
    error('mdemat:analysis:InvalidRefinement', ...
        'Use at least two decreasing positive scales and one error row per scale.');
end
defined = errors(1:end-1, :) > 0 & errors(2:end, :) > 0;
rates = NaN(size(defined));
% Nearby values need log1p; separate logs lose their difference.
logRatios = logRatio(errors(1:end-1, :), errors(2:end, :));
denominators = logRatio(scales(1:end-1), scales(2:end));
candidates = logRatios ./ denominators;
rates(defined) = candidates(defined);
report = struct('Scales', scales, 'Errors', errors, ...
    'LocalRates', rates, 'Defined', defined, 'Finite', isfinite(rates), ...
    'LogScaleRatios', denominators);
end

function result = logRatio(a, b)
result = log(a)-log(b);
near = a > 0 & b > 0 & a/2 <= b & b/2 <= a;
result(near) = log1p((a(near)-b(near))./b(near));
end
