classdef FiniteElementRichardsonMetric < physics.analysis.Metric
    % FINITEELEMENTRICHARDSONMETRIC Richardson extrapolation error metric for
    % finite element spaces.
    %
    %   FiniteElementRichardsonMetric computes relative errors between
    %   solutions on different grid resolutions using Richardson
    %   extrapolation. This approach is invaluable when exact analytical
    %   solutions are unavailable, providing a method to estimate
    %   discretization errors by comparing fine and coarse grid solutions.
    %
    %   The metric implements the Richardson extrapolation principle:
    %   assuming that discretization errors behave as \f$E ≈ Ch^p\f$ where
    %   \f$h\f$ is the grid spacing and \f$p\f$ is the order of
    %   convergence, the difference between fine and coarse grid solutions
    %   provides an estimate of the error.
    %
    %   Richardson extrapolation assumes that the coarse grid is obtained
    %   by uniform coarsening of the fine grid, typically by a factor of 2
    %   in each spatial dimension.
    %
    % See also:
    %   physics.analysis.Metric,
    %   physics.analysis.FiniteElementSpaceAbsoluteMetric,
    %   approx.space.FiniteElementSpace
    
    properties
        Coarse      % Coarse grid finite element space
        Fine        % Fine grid finite element space
        Integrator  % Richardson integrator for nested grid evaluation
        NestEI % Reordered fine element indices for nesting
        NestLI % Reordered fine DOF indices for nesting
    end
    
    properties (Dependent)
        Nodes   % Cell array of integration nodes {fine, coarse}
        Weights % Cell array of integration weights {fine, coarse}
    end
    
    methods
        function obj = FiniteElementRichardsonMetric(fine, coarse, options)
            % FINITEELEMENTRICHARDSONMETRIC Construct an instance of
            % FiniteElementRichardsonMetric.
            %
            %   obj = FiniteElementRichardsonMetric(fine, coarse) creates a
            %   Richardson extrapolation metric for comparing solutions
            %   between @a fine and @a coarse structured grids with empty
            %   reductions list.
            %
            %   obj = FiniteElementRichardsonMetric(fine, coarse,
            %   reductions=reductions) creates a metric with the specified
            %   @a reductions cell array of error reduction types.
            %
            %   The metric requires nested grid structures where the coarse
            %   grid is obtained by uniform coarsening of the fine grid,
            %   enabling Richardson extrapolation error estimates.
            
            arguments
                fine approx.space.FiniteElementSpace
                coarse approx.space.FiniteElementSpace
                options.reductions cell = {}
            end
            
            obj@physics.analysis.Metric(reductions=options.reductions);
            obj.Coarse = coarse;
            obj.Fine = fine;
            
            obj.Integrator = approx.integrate.MultiLevelIntegrator( ...
                coarse.Element.Volume.Rule, nLevels = 2);
            
            obj.setNestElements();
            obj.setNestDofs();
        end
        
        function x = get.Nodes(obj)
            % GET.NODES Get integration nodes for both grid levels.
            
            x = cell(1, 2);
            x{1} = obj.Integrator.Integrators(1).Nodes;
            x{2} = obj.Integrator.Integrators(2).Nodes;
        end
        
        function x = get.Weights(obj)
            % GET.WEIGHTS Get integration weights for both grid levels.
            
            x = cell(1, 2);
            x{1} = obj.Integrator.Integrators(1).Weights;
            x{2} = obj.Integrator.Integrators(2).Weights;
        end
        
        function obj = addReduction(obj, reduction)
            % ADDREDUCTION Add error reduction type to the metric.
            %
            %   obj = addReduction(obj, reduction) adds the specified
            %   error reduction type @a reduction to the metric's list of
            %   reductions to compute during error evaluation.
            
            arguments
                obj physics.analysis.FiniteElementRichardsonMetric
                reduction {mustBeMember(reduction, {'', 'L1', 'L2', 'max'})}
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
                obj physics.analysis.FiniteElementRichardsonMetric
                density {mustBePositive, mustBeInteger}
            end
            
            space = obj.Coarse;
            G = space.Element.Geometry;
            if isa(G, 'core.geometry.Orthotope')
                nd = space.NDims;
                np = ceil(density^(1/nd));
                np = repmat(np, 1, nd);
                obj.Integrator.setPoints(G, np);
            else
                core.except.assert(0, 'InvalidGeometry', ...
                    'FiniteElementRichardsonMetric only supports orthotope geometries.');
            end
        end
        
        function E = eval(obj, U, V)
            % EVAL Compute Richardson extrapolation errors.
            %
            %   E = eval(obj, U, V) computes Richardson extrapolation
            %   errors between fine grid solution @a U and coarse grid
            %   solution @a V using nested grid structure for error
            %   estimation.
            
            arguments
                obj physics.analysis.FiniteElementRichardsonMetric
                U {mustBeNumeric}
                V {mustBeNumeric}
            end
            
            coarseIntegrator = obj.Integrator.Integrators(2);
            nq = length(coarseIntegrator.Weights);
            U = reshape(U, nq, []);
            V = reshape(V, nq, []);
            R = abs(U - V);
            
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
            % REDUCE Apply Richardson extrapolation error reduction.
            
            n = obj.Coarse.NMeshElements;
            h = obj.Coarse.Mesh.computeElementJacobianDeterminants();
            w = obj.Integrator.Integrators(2).Weights;
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
        
        function obj = setNestElements(obj)
            % SETNESTELEMENTS Reorder fine elements to align with coarse
            % elements.
            
            Th0 = obj.Coarse.Mesh;
            Th1 = obj.Fine.Mesh;
            n0 = Th0.Resolution;
            n1 = Th1.Resolution;
            m0 = n1 ./ n0;
            T1 = core.linalg.MultiIndexer(shape=[prod(m0), prod(n0)]);
            I = (1:prod(n1)).';
            LK = T1.linearToMulti(I);
            T2 = core.linalg.MultiIndexer(shape=m0);
            L = T2.linearToMulti(LK(:, 1));
            T3 = core.linalg.MultiIndexer(shape=n0);
            K = T3.linearToMulti(LK(:, 2));
            I = (K-1) .* m0 + (L-1) + 1;
            T4 = core.linalg.MultiIndexer(shape=n1);
            obj.NestEI = T4.multiToLinear(I);
        end
        
        function obj = setNestDofs(obj)
            % SETNESTDOFS Reorder fine DOFs to align with nested structure.
            
            I = obj.NestEI;
            p = obj.Fine.Element.NDofs;
            T = core.linalg.MultiIndexer(shape=[p, size(I, 1)]);
            M = T.factorToMulti({(1:p).', I});
            obj.NestLI = T.multiToLinear(M);
        end
    end
end

