function [F, report] = compress(U, S, V, options)
%COMPRESS Weighted QR/small-core SVD compression without dense expansion.
%   Error is the discarded norm of the represented matrix, not its PDE error.
arguments
    U {mustBeA(U, 'double'), mustBeFinite}
    S {mustBeA(S, 'double'), mustBeFinite}
    V {mustBeA(V, 'double'), mustBeFinite}
    options.LeftMetric = []
    options.RightMetric = []
    options.Axes (1, 2) string = ["x", "v"]
    options.RelativeTolerance (1, 1) double {mustBeNonnegative, mustBeFinite} = 1e-10
    options.AbsoluteTolerance (1, 1) double {mustBeNonnegative, mustBeFinite} = 0
    options.MaxRank (1, 1) double {mustBeNonnegative} = Inf
    options.MinRank (1, 1) double {mustBeNonnegative, mustBeInteger} = 0
end
if ~ismatrix(S) || size(U, 2) ~= size(S, 1) || size(V, 2) ~= size(S, 2)
    error('mdemat:lowrank:Shape', 'Inner factor dimensions must match the core.');
end
if isfinite(options.MaxRank) && fix(options.MaxRank) ~= options.MaxRank
    error('mdemat:lowrank:RankBounds', 'MaxRank must be an integer or Inf.');
end
limit = min([size(U, 1), size(V, 1), options.MaxRank]);
if isnan(limit) || fix(limit) ~= limit || options.MinRank > limit
    error('mdemat:lowrank:RankBounds', 'Rank limits must be ordered nonnegative integers.');
end
[Qx, Rx] = mdemat.common.lowrank.weightedQr(U, options.LeftMetric);
[Qv, Rv] = mdemat.common.lowrank.weightedQr(V, options.RightMetric);
[W, D, Z] = svd(full(Rx * S * Rv'), 'econ'); % Only the small core is dense.
sigma = diag(D);
if options.MinRank > numel(sigma)
    error('mdemat:lowrank:RankBounds', 'MinRank exceeds the available factor subspaces.');
end
tails = sqrt([flipud(cumsum(flipud(abs(sigma).^2))); 0]);
normF = norm(sigma);
threshold = max(options.AbsoluteTolerance, options.RelativeTolerance * normF);
k = find(tails <= threshold, 1) - 1;
if isempty(k), k = numel(sigma); end
k = min([max(k, options.MinRank), limit, numel(sigma)]);
F = mdemat.data.tensor.MatrixFactors(Qx * W(:, 1:k), D(1:k, 1:k), ...
    Qv * Z(:, 1:k), LeftMetric=options.LeftMetric, RightMetric=options.RightMetric, ...
    Axes=options.Axes);
discarded = tails(k + 1);
relative = 0;
if normF > 0, relative = discarded / normF; end
report = struct('Rank', k, 'InputNorm', normF, 'DiscardedNorm', discarded, ...
    'RelativeDiscardedNorm', relative, 'ToleranceMet', discarded <= threshold, ...
    'RankLimitReached', k == limit && discarded > threshold);
end
