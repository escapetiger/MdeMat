function [solution, report] = solveLowRankTransport(initial, velocities, weights, interval, step, options)
%SOLVELOWRANKTRANSPORT Periodic slab transport with isotropic scattering.
%   F(i,j) is a cell average at spatial cell i and velocity j. The equation is
%   F_t + v F_x = Scattering*(<F>-F), <F>=F*w/sum(w).
%   Space uses conservative first-order upwind differences. Time uses either
%   fixed-rank projector splitting or adaptive Euler/compression. The latter
%   requires step*(max(abs(v))/dx+Scattering)<=1 before compression; positivity
%   of a compressed representation is not guaranteed. See the mathematics
%   documentation for the distinction between time, space and rank errors.
arguments
    initial mdemat.data.tensor.MatrixFactors
    velocities (1, :) double {mustBeFinite, mustBeReal}
    weights (:, 1) double {mustBePositive, mustBeFinite, mustBeReal}
    interval (1, 2) double {mustBeFinite}
    step (1, 1) double {mustBePositive, mustBeFinite}
    options.Domain (1, 2) double {mustBeFinite} = [0, 1]
    options.Scattering (1, 1) double {mustBeNonnegative, mustBeFinite} = 0
    options.Method (1, 1) string {mustBeMember(options.Method, ["fixed", "adaptive"])} = "adaptive"
    options.RelativeTolerance (1, 1) double {mustBeNonnegative, mustBeFinite} = 1e-8
    options.MaxRank (1, 1) double {mustBePositive, mustBeInteger} = min(initial.Shape)
    options.PreserveMass (1, 1) logical = true
end
shape = initial.Shape;
nx = shape(1); nv = shape(2);
if nx < 3 || numel(velocities) ~= nv || numel(weights) ~= nv || ...
        options.Domain(2) <= options.Domain(1) || interval(2) < interval(1)
    error('mdemat:radiation:TransportShape', 'Require matching velocity axes, ordered domain/time and at least three cells.');
end
h = diff(options.Domain)/nx;
Mx = h*speye(nx); Mv = spdiags(weights, 0, nv, nv);
if ~isequal(initial.Axes, ["x", "v"]) || ...
        ~isequal(initial.LeftMetric, Mx) || ~isequal(initial.RightMetric, Mv)
    error('mdemat:radiation:TransportMetric', 'Factors must declare cell-volume and velocity quadrature metrics.');
end
if options.Method == "adaptive" && step*(max(abs(velocities))/h+options.Scattering) > 1+10*eps
    error('mdemat:radiation:ExplicitStep', 'The explicit transport/scattering step exceeds its monotonicity bound.');
end
rows = (1:nx)';
back = (speye(nx)-sparse(rows, mod(rows-2, nx)+1, 1, nx, nx))/h;
forward = (sparse(rows, mod(rows, nx)+1, 1, nx, nx)-speye(nx))/h;
left = {-back, -forward, options.Scattering*speye(nx)};
right = {diag(max(velocities, 0)), diag(min(velocities, 0)), ...
    ones(nv, 1)*weights'/sum(weights)-eye(nv)};
rhs = @(~, F) mdemat.common.operator.applySeparated(left, right, F);
started = tic;
if options.Method == "fixed"
    [field, evolution] = mdemat.common.ode.evolveFixedRank(initial, rhs, interval, step);
else
    [field, evolution] = mdemat.common.ode.evolveAdaptive(initial, rhs, interval, step, ...
        RelativeTolerance=options.RelativeTolerance, MaxRank=options.MaxRank, ...
        PreserveMass=options.PreserveMass);
end
duration = toc(started);
solution = struct('Field', field, 'Time', interval(2), ...
    'CellCenters', options.Domain(1)+((1:nx)-0.5)*h, ...
    'Velocities', velocities, 'VelocityWeights', weights);
report = struct('Status', 'COMPLETED', 'Assessment', 'UNVERIFIED', ...
    'SpatialOrder', 1, 'Evolution', evolution, 'Duration', duration, ...
    'InitialMass', mass(initial, weights, h), 'FinalMass', mass(field, weights, h), ...
    'DenseEntries', prod(shape), 'FactorEntries', field.StoredEntries, ...
    'PositivityGuaranteedAfterCompression', false);
end

function value = mass(F, w, h)
value = h*(sum(F.U, 1)*F.S*(F.V'*w));
end
