classdef MatrixFactors
    %MATRIXFACTORS Optional matrix representation F=U*S*V' with explicit metrics.
    %   Native dense/sparse paths do not require this value object. materialize
    %   is explicit; sampling and norms do not expand the complete field.
    properties (SetAccess = private)
        U
        S
        V
        LeftMetric
        RightMetric
        Axes
    end
    properties (Dependent)
        Shape
        Rank
        StoredEntries
    end
    methods
        function obj = MatrixFactors(U, S, V, options)
            arguments
                U {mustBeA(U, 'double'), mustBeFinite}
                S {mustBeA(S, 'double'), mustBeFinite}
                V {mustBeA(V, 'double'), mustBeFinite}
                options.LeftMetric = []
                options.RightMetric = []
                options.Axes (1, 2) string = ["x", "v"]
            end
            r = size(S, 1);
            if ~ismatrix(U) || ~ismatrix(S) || ~ismatrix(V) || ...
                    size(S, 2) ~= r || size(U, 2) ~= r || size(V, 2) ~= r
                error('mdemat:tensor:FactorShape', 'Factors must have shapes (nx,r), (r,r), (nv,r).');
            end
            validateMetric(options.LeftMetric, size(U, 1));
            validateMetric(options.RightMetric, size(V, 1));
            obj.U = U; obj.S = S; obj.V = V;
            obj.LeftMetric = options.LeftMetric;
            obj.RightMetric = options.RightMetric;
            obj.Axes = options.Axes;
        end
        function value = get.Shape(obj)
            value = [size(obj.U, 1), size(obj.V, 1)];
        end
        function value = get.Rank(obj)
            % Representation width, not a hidden numerical-rank calculation.
            value = size(obj.S, 1);
        end
        function value = get.StoredEntries(obj)
            % Factor payload only; metrics/operators/workspaces are separate.
            value = numel(obj.U) + numel(obj.S) + numel(obj.V);
        end
        function F = materialize(obj)
            F = obj.U * obj.S * obj.V';
        end
        function result = withFactors(obj, U, S, V)
            % Reuse immutable, already-validated metrics on the same axes.
            validateattributes(U, {'double'}, {'2d', 'finite'});
            validateattributes(S, {'double'}, {'2d', 'finite'});
            validateattributes(V, {'double'}, {'2d', 'finite'});
            r = size(S, 1);
            if size(S, 2) ~= r || size(U, 2) ~= r || size(V, 2) ~= r || ...
                    ~isequal([size(U, 1), size(V, 1)], obj.Shape)
                error('mdemat:tensor:FactorShape', 'Replacement factors must preserve field axes.');
            end
            result = obj;
            result.U = U; result.S = S; result.V = V;
        end
        function values = sample(obj, rows, columns)
            values = obj.U(rows, :) * obj.S * obj.V(columns, :)';
        end
        function value = weightedNorm(obj)
            Gx = obj.U' * metricAction(obj.LeftMetric, obj.U);
            Gv = obj.V' * metricAction(obj.RightMetric, obj.V);
            value = sqrt(max(0, real(trace(obj.S' * Gx * obj.S * Gv))));
        end
        function record = snapshot(obj)
            % Portable completed-step representation; no live object graph.
            record = struct('Schema', 'mdemat.matrix-factors', 'Version', 1, ...
                'U', obj.U, 'S', obj.S, 'V', obj.V, ...
                'LeftMetric', obj.LeftMetric, 'RightMetric', obj.RightMetric, ...
                'Axes', obj.Axes);
        end
    end
    methods (Static)
        function obj = fromSnapshot(record)
            required = {'Schema', 'Version', 'U', 'S', 'V', ...
                'LeftMetric', 'RightMetric', 'Axes'};
            if ~isstruct(record) || ~isscalar(record) || ...
                    ~isequal(sort(fieldnames(record)), sort(required(:))) || ...
                    ~isequal(record.Schema, 'mdemat.matrix-factors') || ...
                    ~isequal(record.Version, 1)
                error('mdemat:tensor:SnapshotSchema', 'Unsupported matrix-factor snapshot.');
            end
            obj = mdemat.data.tensor.MatrixFactors(record.U, record.S, record.V, ...
                LeftMetric=record.LeftMetric, RightMetric=record.RightMetric, ...
                Axes=record.Axes);
        end
    end
end

function validateMetric(M, n)
if isempty(M)
    return;
end
validateattributes(M, {'double'}, {'2d', 'finite', 'size', [n, n]});
if ~ishermitian(M)
    error('mdemat:tensor:Metric', 'Metric must be Hermitian positive definite.');
end
[~, flag] = chol(M);
if flag ~= 0
    error('mdemat:tensor:Metric', 'Metric must be positive definite.');
end
end

function value = metricAction(M, X)
if isempty(M), value = X; else, value = M * X; end
end
