function report = weightedError(values, reference, weights, options)
%WEIGHTEDERROR Physical norms of finite components-by-points samples.
%   R = weightedError(U,V,W) uses nonnegative physical quadrature weights,
%   including Jacobians. Zero-weight points are excluded. No volume or
%   component normalization is applied. Complex samples use their magnitude.
%   ComponentLp contains one spatial norm per component. L1 sums ComponentL1,
%   L2 is the Euclidean norm of ComponentL2, and LInf is max(ComponentLInf).
%   Callers form spatial-then-component Euclidean norms explicitly.
%
%   Samples are promoted to double. Logarithmic scaling prevents overflowing
%   intermediate differences from destroying representable weighted norms.
%   Unrepresentable output norms remain Inf. RelativeL2 is NaN unless both L2
%   norms and their ratio are finite and ReferenceL2 > referenceFloor.
arguments
    values {mustBeFloat, mustBeFinite}
    reference {mustBeFloat, mustBeFinite}
    weights {mustBeNumeric, mustBeReal, mustBeVector, mustBeFinite, mustBeNonnegative}
    options.referenceFloor (1, 1) double {mustBeNonnegative, mustBeFinite} = 0
end
if ~ismatrix(values) || ~isequal(size(values), size(reference)) || ...
        size(values, 1) == 0 || size(values, 2) ~= numel(weights) || ~any(weights > 0)
    error('mdemat:analysis:InvalidShape', ...
        'Values/reference must be components-by-points with positive total weight.');
end
active = weights > 0;
weights = double(weights(active));
weights = weights(:).';
values = double(values(:, active));
reference = double(reference(:, active));
logError = logMagnitudeDifference(values, reference);
logReference = logMagnitude(reference);
logWeights = log(weights);
componentLogL1 = logSum(logError + logWeights, 2);
componentLogL2 = logSum(2*logError + logWeights, 2)/2;
componentLogLInf = max(logError, [], 2);
logL1 = logSum(componentLogL1, 1);
logL2 = logSum(2*componentLogL2, 1)/2;
logReferenceL2 = logSum(logSum(2*logReference + logWeights, 2), 1)/2;
l1 = exp(logL1);
l2 = exp(logL2);
linf = exp(max(componentLogLInf));
referenceL2 = exp(logReferenceL2);
% Keep ordinary arithmetic accuracy when its intermediates are representable.
% Log scaling is the fallback, not a reason to perturb small exact integers.
delta = values-reference;
weighted = delta.*sqrt(weights);
absoluteDelta = abs(delta);
componentL1 = exp(componentLogL1);
componentL2 = exp(componentLogL2);
componentLInf = exp(componentLogLInf);
if all(isfinite(absoluteDelta), 'all') && all(isfinite(weighted), 'all')
    directL1 = absoluteDelta*weights.';
    directL2 = arrayfun(@(row) norm(weighted(row,:), 2), (1:size(values,1)).');
    directLInf = max(absoluteDelta, [], 2);
    use = directL1 > 0 | isinf(componentLogL1);
    componentL1(use) = directL1(use);
    use = directL2 > 0 | isinf(componentLogL2);
    componentL2(use) = directL2(use);
    componentLInf = directLInf;
    l1 = sum(componentL1);
    l2 = norm(componentL2, 2);
    linf = max(componentLInf);
end
weightedReference = reference.*sqrt(weights);
if all(isfinite(weightedReference), 'all')
    directReference = norm(weightedReference, 'fro');
    if directReference > 0 || isinf(logReferenceL2), referenceL2 = directReference; end
end
relative = NaN;
defined = isfinite(l2) && isfinite(referenceL2) && referenceL2 > options.referenceFloor;
if defined
    relative = l2/referenceL2;
    defined = isfinite(relative);
    if ~defined, relative = NaN; end
end
report = struct('L1', l1, 'L2', l2, 'LInf', linf, ...
    'ReferenceL2', referenceL2, 'RelativeL2', relative, ...
    'RelativeDefined', defined, 'RelativeFinite', isfinite(relative), ...
    'ReferenceFloor', options.referenceFloor, ...
    'ComponentL1', componentL1, 'ComponentL2', componentL2, ...
    'ComponentLInf', componentLInf, ...
    'FiniteNorms', all(isfinite([l1, l2, linf, referenceL2])), ...
    'Norm', 'physical weighted samples; explicit component aggregation');
end

function result = logMagnitudeDifference(a, b)
result = combine(logRealDifference(real(a), real(b)), ...
    logRealDifference(imag(a), imag(b)));
end

function result = logRealDifference(a, b)
difference = a-b;
result = log(abs(difference));
overflow = isinf(difference);
if any(overflow, 'all')
    scale = max(abs(a(overflow)), abs(b(overflow)));
    result(overflow) = log(scale) + log(abs(a(overflow)./scale-b(overflow)./scale));
end
end

function result = logMagnitude(a)
result = combine(log(abs(real(a))), log(abs(imag(a))));
end

function result = combine(a, b)
anchor = max(a, b);
zero = isinf(anchor) & anchor < 0;
anchor(zero) = 0;
result = anchor + log(exp(2*(a-anchor))+exp(2*(b-anchor)))/2;
result(zero) = -Inf;
end

function result = logSum(values, dimension)
anchor = max(values, [], dimension);
anchor(isinf(anchor) & anchor < 0) = 0;
result = anchor + log(sum(exp(values-anchor), dimension));
end
