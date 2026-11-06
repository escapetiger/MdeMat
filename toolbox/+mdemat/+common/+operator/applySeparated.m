function [result, report] = applySeparated(left, right, field, options)
%APPLYSEPARATED Apply sum_j left{j}*F*right{j}' to a factorized field.
%   Contractions are formed from factors and a block core; the full field and
%   global Kronecker operator are not materialized. Compression error is explicit.
arguments
    left cell
    right cell
    field mdemat.data.tensor.MatrixFactors
    options.RelativeTolerance (1, 1) double {mustBeNonnegative} = 0
    options.MaxRank (1, 1) double {mustBeNonnegative} = Inf
end
if numel(left) ~= numel(right)
    error('mdemat:operator:Terms', 'Left/right operator term counts must agree.');
end
shape = field.Shape;
U = zeros(shape(1), 0); V = zeros(shape(2), 0); S = zeros(0);
for term = 1:numel(left)
    validateattributes(left{term}, {'double'}, {'finite', 'size', [shape(1), shape(1)]});
    validateattributes(right{term}, {'double'}, {'finite', 'size', [shape(2), shape(2)]});
    U = [U, left{term} * field.U]; %#ok<AGROW>
    V = [V, right{term} * field.V]; %#ok<AGROW>
    S = blkdiag(S, field.S);
end
[result, report] = mdemat.common.lowrank.compress(U, S, V, ...
    LeftMetric=field.LeftMetric, RightMetric=field.RightMetric, ...
    Axes=field.Axes, RelativeTolerance=options.RelativeTolerance, MaxRank=options.MaxRank);
end
