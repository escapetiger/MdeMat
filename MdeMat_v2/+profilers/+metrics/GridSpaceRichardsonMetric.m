classdef GridSpaceRichardsonMetric < profilers.metrics.Metric
    % GRIDSPACERICHARDSONMETRIC Richardson extrapolation error metric for
    % grid spaces.
    %
    %   GridSpaceRichardsonMetric computes relative errors between
    %   solutions on different grid resolutions using Richardson
    %   extrapolation. This approach is invaluable when exact analytical
    %   solutions are unavailable, providing a method to estimate
    %   discretization errors by comparing fine and coarse grid solutions.
    %
    %   The metric implements the Richardson extrapolation principle:
    %   assuming that discretization errors behave as E ≈ Ch^p where h is
    %   the grid spacing and p is the order of convergence, the difference
    %   between fine and coarse grid solutions provides an estimate of the
    %   error.
    %
    % Examples:
    %   % Create Richardson metric for L2 error estimation
    %   metric = GridSpaceRichardsonMetric({'L2'}, fineSpace, coarseSpace);
    %   metric.setPoints([6, 6]); % 6×6 integration points per element
    %   
    %   % Compute Richardson extrapolation error
    %   fineNodes = metric.nodes{1};
    %   coarseNodes = metric.nodes{2};
    %   fineValues = fineSpace.evaluate(fineNodes, fineCoeffs);
    %   coarseValues = coarseSpace.evaluate(coarseNodes, coarseCoeffs);
    %   error = metric.evaluate(fineValues, coarseValues);
    %
    %   % Use nested DOF ordering for convenience
    %   reorderedCoeffs = fineCoeffs(metric.iNestDofs, :);
    %   error = metric.evaluate(fineValues, coarseValues);
    %
    % Notes:
    %   Richardson extrapolation assumes that the coarse grid is obtained
    %   by uniform coarsening of the fine grid, typically by a factor of 2
    %   in each spatial dimension.
    %
    % See also:
    %   profilers.metrics.Metric, profilers.metrics.MeshSpaceAbsoluteMetric,
    %   approx.space.MeshSpace
    
    properties
        coarse      % Coarse grid finite element space
        fine        % Fine grid finite element space
        integrator  % Richardson integrator for nested grid evaluation

        iNestElements % Reordered fine element indices for nesting
        iNestDofs     % Reordered fine DOF indices for nesting
    end

    properties (Dependent)
        nodes   % Cell array of integration nodes {fine, coarse}
        weights % Cell array of integration weights {fine, coarse}
    end
    
    methods
        function obj = GridSpaceRichardsonMetric(reduction, fine, coarse)
            % GRIDSPACERICHARDSONMETRIC Constructor for Richardson metric.
            %
            %   obj = GridSpaceRichardsonMetric(reduction, fine, coarse)
            %   creates a Richardson extrapolation metric for comparing
            %   solutions between fine and coarse structured grids.
            %
            % Inputs:
            %   reduction - Cell array of error reduction types
            %   fine - Fine grid MeshSpace object
            %   coarse - Coarse grid MeshSpace object (must be nested)
            %
            % Outputs:
            %   obj - Constructed GridSpaceRichardsonMetric object

            core.except.assert(isa(fine, 'approx.space.MeshSpace'), ...
                'InvalidInput', 'Fine space must be a mesh space.');
            core.except.assert(isa(coarse, 'approx.space.MeshSpace'), ...
                'InvalidInput', 'Coarse space must be a mesh space.');
            core.except.assert(isa(fine.mesh, 'approx.mesh.Grid'), ...
                'InvalidInput', 'Fine mesh must be a grid.');
            core.except.assert(isa(coarse.mesh, 'approx.mesh.Grid'), ...
                'InvalidInput', 'Coarse mesh must be a grid.');

            obj@profilers.metrics.Metric(reduction);
            obj.coarse = coarse;
            obj.fine = fine;
            rule = coarse.element.volumeData.integrator.rule;
            obj.integrator = approx.integrate.RichardsonIntegrator(rule, 2);
            obj.setNestElements();
            obj.setNestDofs();
        end

        function x = get.nodes(obj)
            % GET.NODES Get integration nodes for both grid levels.
            %
            %   x = get.nodes(obj) returns a cell array containing
            %   integration nodes for fine and coarse grid evaluation.
            %
            % Inputs:
            %   obj - The GridSpaceRichardsonMetric object
            %
            % Outputs:
            %   x - {fine_nodes, coarse_nodes} cell array

            x = obj.integrator.nodes;
        end

        function x = get.weights(obj)
            % GET.WEIGHTS Get integration weights for both grid levels.
            %
            %   x = get.weights(obj) returns a cell array containing
            %   integration weights for fine and coarse grid evaluation.
            %
            % Inputs:
            %   obj - The GridSpaceRichardsonMetric object
            %
            % Outputs:
            %   x - {fine_weights, coarse_weights} cell array

            x = obj.integrator.weights;
        end

        function obj = setPoints(obj, n)
            % SETPOINTS Configure integration rules for both grid levels.
            %
            %   obj = setPoints(obj, n) sets the integration rule to use n
            %   quadrature points per spatial dimension for both fine and
            %   coarse grid evaluations in the Richardson extrapolation.
            %
            % Inputs:
            %   obj - The GridSpaceRichardsonMetric object
            %   n - Number or vector of integration points per dimension
            %
            % Outputs:
            %   obj - The GridSpaceRichardsonMetric object
            
            G = obj.coarse.element.geometry;
            obj.integrator.setPoints(n, G.lower, G.upper);
        end

        function E = evaluate(obj, U, V)
            % EVALUATE Compute Richardson extrapolation errors.
            %
            %   E = evaluate(obj, U, V) computes Richardson extrapolation
            %   errors between fine grid solution U and coarse grid
            %   solution V. The method accounts for proper domain
            %   integration and element scaling to provide meaningful error
            %   estimates.
            %
            % Inputs:
            %   obj - The GridSpaceRichardsonMetric object
            %   U - Fine grid solution values at integration points
            %   V - Coarse grid solution values at integration points
            %
            % Outputs:
            %   E - Matrix of computed Richardson errors for each reduction type

            I = obj.integrator;
            w = I.weights{2};
            U = reshape(U, length(w), []);
            V = reshape(V, length(w), []);

            %< Compute absolute difference
            R = abs(U - V);
            
            %< Get mesh properties for proper scaling
            Th = obj.coarse.mesh;
            n = Th.nTotalElements;
            h = Th.magnitudes;

            %< Compute each requested error norm
            nTypes = length(obj.reduction);
            E = cell(nTypes, 1);
            for i = 1:nTypes
                red = obj.reduction{i};
                switch red
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
                    case 'Lx'
                        e = max(R, [], 1);
                        e = reshape(e, n, []);
                        e = max(e, [], 1);
                    case 'none'
                        e = R;
                end
                E{i} = e(:).';
            end
            E = cell2mat(E);
        end
    end

    methods (Access = protected)
        function obj = setNestElements(obj)
            % SETNESTELEMENTS Reorder fine elements to align with coarse
            % elements.
            %
            %   obj = setNestElements(obj) computes the reordering of fine
            %   grid elements to align with the nested coarse grid
            %   structure. This enables proper Richardson extrapolation
            %   between grids.
            %
            % Inputs:
            %   obj - The GridSpaceRichardsonMetric object
            %
            % Outputs:
            %   obj - The GridSpaceRichardsonMetric object

            Th0 = obj.coarse.mesh;
            Th1 = obj.fine.mesh;
            n0 = Th0.resolution;
            n1 = Th1.resolution;
            m0 = n1 ./ n0;
            T1 = core.linalg.MultiIndexer([prod(m0), prod(n0)]);
            I = (1:prod(n1)).';
            LK = T1.linearToMulti(I);
            T2 = core.linalg.MultiIndexer(m0);
            L = T2.linearToMulti(LK(:, 1));
            T3 = core.linalg.MultiIndexer(n0);
            K = T3.linearToMulti(LK(:, 2));
            I = (K-1) .* m0 + (L-1) + 1;
            T4 = core.linalg.MultiIndexer(n1);
            obj.iNestElements = T4.multiToLinear(I);
        end

        function obj = setNestDofs(obj)
            % SETNESDDOFS Reorder fine DOFs to align with nested structure.
            %
            %   obj = setNestDofs(obj) computes the reordering of fine grid
            %   degrees of freedom to align with the nested element
            %   structure, facilitating coefficient reordering.
            %
            % Inputs:
            %   obj - The GridSpaceRichardsonMetric object
            %
            % Outputs:
            %   obj - The GridSpaceRichardsonMetric object

            I = obj.iNestElements;
            p = obj.fine.element.nDofs;
            T = core.linalg.MultiIndexer([p, size(I, 1)]);
            M = T.factorToMulti({(1:p).', I});
            obj.iNestDofs = T.multiToLinear(M);
        end
    end
end