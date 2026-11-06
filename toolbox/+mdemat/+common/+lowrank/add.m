function [F, report] = add(A, B, options)
%ADD Sum two represented fields and recompress with an explicit error budget.
arguments
    A mdemat.data.tensor.MatrixFactors
    B mdemat.data.tensor.MatrixFactors
    options.ScaleB (1, 1) double {mustBeFinite} = 1
    options.RelativeTolerance (1, 1) double {mustBeNonnegative} = 0
    options.AbsoluteTolerance (1, 1) double {mustBeNonnegative} = 0
    options.MaxRank (1, 1) double {mustBeNonnegative} = Inf
end
if ~isequal(A.Shape, B.Shape) || ~isequal(A.Axes, B.Axes) || ...
        ~isequal(A.LeftMetric, B.LeftMetric) || ~isequal(A.RightMetric, B.RightMetric)
    error('mdemat:lowrank:IncompatibleFields', 'Field axes and metrics must agree.');
end
[F, report] = mdemat.common.lowrank.compress([A.U, B.U], ...
    blkdiag(A.S, options.ScaleB * B.S), [A.V, B.V], ...
    LeftMetric=A.LeftMetric, RightMetric=A.RightMetric, ...
    Axes=A.Axes, ...
    RelativeTolerance=options.RelativeTolerance, AbsoluteTolerance=options.AbsoluteTolerance, ...
    MaxRank=options.MaxRank);
end
