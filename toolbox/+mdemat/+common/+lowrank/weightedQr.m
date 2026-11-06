function [Q, R] = weightedQr(A, metric)
%WEIGHTEDQR Factor A=Q*R with Q'*metric*Q=I for an SPD physical metric.
%   Empty metric means the Euclidean inner product. Dependent columns are
%   represented by a rank-deficient R; no inverse of a small singular core is used.
arguments
    A {mustBeA(A, 'double'), mustBeFinite}
    metric = []
end
if ~ismatrix(A)
    error('mdemat:lowrank:Shape', 'Expected a matrix.');
end
if isempty(metric)
    [Q, R] = qr(A, 0);
else
    validateattributes(metric, {'double'}, {'finite', 'size', [size(A, 1), size(A, 1)]});
    if ~ishermitian(metric)
        error('mdemat:lowrank:Metric', 'Metric must be Hermitian positive definite.');
    end
    [C, flag] = chol(metric);
    if flag ~= 0
        error('mdemat:lowrank:Metric', 'Metric must be positive definite.');
    end
    [W, R] = qr(C * A, 0);
    Q = C \ W;
end
end
