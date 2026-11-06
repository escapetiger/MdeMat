function [field, report] = evolveAdaptive(initial, rhs, interval, step, options)
%EVOLVEADAPTIVE Rank-adaptive explicit Euler followed by weighted compression.
%   This is a retraction/compression method, not a projector-splitting method.
%   PreserveMass corrects the scalar weighted integral inside retained factor
%   subspaces. It fails explicitly if those spaces cannot represent the moment.
%   Positivity/AP/stiff stability are not inferred from mass preservation.
arguments
    initial mdemat.data.tensor.MatrixFactors
    rhs function_handle
    interval (1, 2) double {mustBeFinite}
    step (1, 1) double {mustBePositive, mustBeFinite}
    options.RelativeTolerance (1, 1) double {mustBeNonnegative} = 1e-8
    options.AbsoluteTolerance (1, 1) double {mustBeNonnegative} = 0
    options.MaxRank (1, 1) double {mustBePositive, mustBeInteger} = min(initial.Shape)
    options.PreserveMass (1, 1) logical = false
end
if interval(2) < interval(1)
    error('mdemat:ode:TimeInterval', 'Time interval must be ordered.');
end
field = initial; t = interval(1);
times = t; ranks = field.Rank; errors = 0; masses = integral(field);
while t < interval(2)
    h = min(step, interval(2)-t);
    derivative = rhs(t, field);
    if ~isa(derivative, 'mdemat.data.tensor.MatrixFactors')
        error('mdemat:ode:FactorRhs', 'RHS must return a factorized coefficient field.');
    end
    [raw, ~] = mdemat.common.lowrank.add(field, derivative, ScaleB=h);
    targetMass = integral(raw);
    [field, compression] = mdemat.common.lowrank.compress(raw.U, raw.S, raw.V, ...
        LeftMetric=raw.LeftMetric, RightMetric=raw.RightMetric, ...
        Axes=raw.Axes, ...
        RelativeTolerance=options.RelativeTolerance, AbsoluteTolerance=options.AbsoluteTolerance, ...
        MaxRank=options.MaxRank);
    correctionNorm = 0;
    if options.PreserveMass
        [field, correctionNorm] = correctMass(field, targetMass);
    end
    t = t + h;
    times(end+1) = t; %#ok<AGROW>
    ranks(end+1) = field.Rank; %#ok<AGROW>
    errors(end+1) = compression.DiscardedNorm + correctionNorm; %#ok<AGROW>
    masses(end+1) = integral(field); %#ok<AGROW>
end
report = struct('Method', 'rank-adaptive-Euler', 'Order', 1, 'Times', times, ...
    'Ranks', ranks, 'RepresentationChangeBounds', errors, 'Masses', masses, ...
    'StoredEntries', field.StoredEntries);
end

function value = integral(F)
left = ones(F.Shape(1), 1); right = ones(F.Shape(2), 1);
value = (left' * action(F.LeftMetric, F.U)) * F.S ...
    * (F.V' * action(F.RightMetric, right));
end

function [F, change] = correctMass(F, target)
left = F.U' * action(F.LeftMetric, ones(F.Shape(1), 1));
right = F.V' * action(F.RightMetric, ones(F.Shape(2), 1));
denominator = norm(left)^2 * norm(right)^2;
delta = target - integral(F);
if abs(delta) <= 100*eps*max(1, abs(target))
    change = 0;
    return;
end
if denominator <= realmin
    error('mdemat:ode:MassSubspaceLost', 'Retained factor spaces cannot preserve the requested mass.');
end
correction = delta * (left * right') / denominator;
F = F.withFactors(F.U, F.S + correction, F.V);
change = norm(correction, 'fro');
end

function value = action(M, X)
if isempty(M), value = X; else, value = M*X; end
end
