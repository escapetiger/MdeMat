classdef FiniteElementAbsoluteMetric < physics.analysis.Metric
    % FINITEELEMENTABSOLUTEMETRIC Absolute error computation in finite
    % element spaces.
    %
    %   FiniteElementAbsoluteMetric computes absolute errors between
    %   numerical finite element solutions and reference solutions (either
    %   exact analytical solutions or other numerical solutions). This
    %   class is fundamental for finite element error analysis and
    %   convergence studies.
    %
    %   The metric integrates errors over the computational domain using
    %   appropriate quadrature rules, accounting for element sizes and
    %   geometric transformations. It supports multiple error norms
    %   simultaneously for comprehensive error analysis.
    %
    % See also:
    %   physics.analysis.Metric, approx.space.FiniteElementSpace,
    %   physics.analysis.FiniteElementRichardsonMetric

    properties
        Space % Finite element mesh space for error computation
        Integrator % Numerical integrator for error integration
    end

    properties (Dependent)
        Nodes % Integration nodes on reference element
        Weights % Integration weights on reference element
    end

    methods
        function obj = FiniteElementAbsoluteMetric(space, options)
            % FINITEELEMENTABSOLUTEMETRIC Construct an instance of
            % FiniteElementAbsoluteMetric.
            %
            %   obj = FiniteElementAbsoluteMetric(space) creates an
            %   absolute error metric for the specified finite element mesh
            %   space @a space with empty reductions list.
            %
            %   obj = FiniteElementAbsoluteMetric(space,
            %   reductions=reductions) creates an absolute error metric
            %   with the specified @a reductions cell array of error
            %   reduction types.

            arguments
                space approx.space.FiniteElementSpace
                options.reductions cell = {}
            end

            obj@physics.analysis.Metric(reductions = options.reductions);
            obj.Space = space;
            obj.Integrator = space.Element.Volume.copy();
        end

        function X = get.Nodes(obj)
            % GET.NODES Get integration nodes on the reference element.

            X = obj.Integrator.Nodes;
        end

        function w = get.Weights(obj)
            % GET.WEIGHTS Get integration weights on the reference element.

            w = obj.Integrator.Weights;
        end

        function obj = addReduction(obj, reduction)
            % ADDREDUCTION Add error reduction type to the metric.
            %
            %   obj = addReduction(obj, reduction) adds the specified
            %   error reduction type @a reduction to the metric's list of
            %   reductions to compute during error evaluation.

            arguments
                obj physics.analysis.FiniteElementAbsoluteMetric
                reduction{mustBeMember(reduction, {'', 'L1', 'L2', 'max'})}
            end

            obj.Reductions{end+1} = reduction;
        end

        function obj = setPoints(obj, density)
            % SETPOINTS Configure integration nodes and weights.
            %
            %   obj = setPoints(obj, density) sets the number of
            %   integration points per element for the quadrature rule used
            %   in error integration.

            arguments
                obj physics.analysis.FiniteElementAbsoluteMetric
                density{mustBePositive, mustBeInteger}
            end

            space = obj.Space;
            G = space.Element.Geometry;
            if isa(G, 'core.geometry.Orthotope')
                nd = space.NDims;
                np = ceil(density^(1 / nd));
                np = repmat(np, 1, nd);
                obj.Integrator.setPoints(G, np);
            else
                core.except.assert(0, 'InvalidGeometry', ...
                    'FiniteElementAbsoluteMetric only supports orthotope geometries.');
            end
        end

        function E = eval(obj, U, V)
            % EVAL Compute absolute errors between solutions.
            %
            %   E = eval(obj, U, V) computes absolute errors between
            %   numerical solution @a U and reference solution @a V using
            %   domain integration with proper element scaling.

            arguments
                obj physics.analysis.FiniteElementAbsoluteMetric
                U {mustBeNumeric}
                V {mustBeNumeric}
            end

            nq = obj.Integrator.NPoints;
            U = reshape(U, nq, []);
            V = reshape(V, nq, []);
            R = abs(V-U);

            nt = length(obj.Reductions);
            E = cell(nt, 1);
            for i = 1:nt
                E{i} = obj.reduce(R, obj.Reductions{i});
            end
            E = cell2mat(E);
        end
    end

    methods (Access = protected)
        function e = reduce(obj, R, reduction)
            % REDUCE Apply error reduction over the computational domain.

            n = obj.Space.NMeshElements;
            h = obj.Space.Mesh.computeElementJacobianDeterminants();
            w = obj.Integrator.Weights;
            switch reduction
                case 'L1'
                    e = w * R;
                    e = reshape(e, n, []);
                    e = e .* h(:);
                    e = sum(e, 1);
                case 'L2'
                    e = w * R.^2;
                    e = reshape(e, n, []);
                    e = e .* h(:);
                    e = sqrt(sum(e, 1));
                case 'max'
                    e = max(R, [], 1);
                    e = reshape(e, n, []);
                    e = max(e, [], 1);
                otherwise
                    e = R;
            end
            e = e(:).';
        end
    end
end

function mustBeFunctionOrEmpty(x)
if ~isempty(x)
    mustBeA(x, 'function_handle');
end
end