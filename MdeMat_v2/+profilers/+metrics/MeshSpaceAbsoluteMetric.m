classdef MeshSpaceAbsoluteMetric < profilers.metrics.Metric
    % MESHSPACEABSOLUTEMETRIC Absolute error computation in finite element
    % spaces.
    %
    %   MeshSpaceAbsoluteMetric computes absolute errors between numerical
    %   finite element solutions and reference solutions (either exact
    %   analytical solutions or other numerical solutions). This class is
    %   fundamental for finite element error analysis and convergence
    %   studies.
    %
    %   The metric integrates errors over the computational domain using
    %   appropriate quadrature rules, accounting for element sizes and
    %   geometric transformations. It supports multiple error norms
    %   simultaneously for comprehensive error analysis.
    %
    % Examples:
    %   % Create metric for L2 error analysis
    %   metric = MeshSpaceAbsoluteMetric({'L2'}, meshSpace);
    %   metric.setPoints([8, 8]); % 8×8 integration points per element
    %   
    %   % Compute error against exact solution
    %   exactFunc = @(x) sin(pi*x(1,:)) .* cos(pi*x(2,:));
    %   nodes = metric.nodes;
    %   exactValues = meshSpace.evaluate(exactFunc, nodes);
    %   numericalValues = meshSpace.evaluate(nodes, coefficients);
    %   error = metric.evaluate(numericalValues, exactValues);
    %
    %   % Multiple error norms for comprehensive analysis
    %   metric = MeshSpaceAbsoluteMetric({'L1', 'L2', 'Lx'}, meshSpace);
    %   errors = metric.evaluate(numerical, exact);
    %
    % See also:
    %   profilers.metrics.Metric, approx.space.MeshSpace,
    %   profilers.metrics.GridSpaceRichardsonMetric

    properties
        space       % Finite element mesh space for error computation
        integrator  % Numerical integrator for error integration
    end

    properties (Dependent)
        nodes   % Integration nodes on reference element
        weights % Integration weights on reference element
    end

    methods
        function obj = MeshSpaceAbsoluteMetric(reduction, space)
            % MESHSPACEABSOLUTEMETRIC Constructor for mesh space absolute
            % metric.
            %
            %   obj = MeshSpaceAbsoluteMetric(reduction, space) creates an
            %   absolute error metric for the specified finite element mesh
            %   space with the given error reduction types.
            %
            % Inputs:
            %   reduction - Cell array of error reduction types
            %   space - MeshSpace object defining the finite element space
            %
            % Outputs:
            %   obj - Constructed MeshSpaceAbsoluteMetric object

            core.except.assert(isa(space, 'approx.space.MeshSpace'), ...
                'InvalidInput', 'Input must be a mesh space.');

            obj@profilers.metrics.Metric(reduction);
            obj.space = space;
            obj.integrator = space.element.volumeData.integrator.copy();
        end

        function X = get.nodes(obj)
            % GET.NODES Get integration nodes on the reference element.

            X = obj.integrator.nodes;
        end

        function w = get.weights(obj)
            % GET.WEIGHTS Get integration weights on the reference element.

            w = obj.integrator.weights;
        end

        function obj = setPoints(obj, n)
            % SETPOINTS Configure integration nodes and weights.
            %
            %   obj = setPoints(obj, n) sets the integration rule to use n
            %   quadrature points per spatial dimension. Higher values
            %   provide more accurate error integration but increase
            %   computational cost.
            %
            % Inputs:
            %   obj - The MeshSpaceAbsoluteMetric object
            %   n - Number or vector of integration points per dimension
            %
            % Outputs:
            %   obj - The MeshSpaceAbsoluteMetric object

            G = obj.space.element.geometry;
            obj.integrator.setPoints(n, G.lower, G.upper);
        end

        function E = evaluate(obj, U, varargin)
            % EVALUATE Compute absolute errors between solutions.
            %
            %   E = evaluate(obj, U, V) computes absolute errors between
            %   numerical solution U and reference solution V using domain
            %   integration with proper element scaling.
            %
            %   E = evaluate(obj, U, f, ...) computes errors against a
            %   function handle f evaluated with additional arguments.
            %
            % Inputs:
            %   obj - The MeshSpaceAbsoluteMetric object
            %   U - Numerical solution values at integration points
            %   varargin - Input arguments
            %<   V - Reference solution values (alternative 1)
            %<   f - Reference solution function handle (alternative 2)
            %
            % Outputs:
            %   E - Matrix of computed errors, one row per reduction type

            I = obj.integrator;
            m = I.nPoints;
            w = I.weights;
            U = reshape(U, m, []);

            %< Evaluate reference solution
            if length(varargin) == 1
                V = reshape(varargin{1}, m, []);
            else
                f = varargin{1};
                V = reshape(f(varargin{2:end}), m, []);
            end

            %< Get mesh properties for proper scaling
            Th = obj.space.mesh;
            n = Th.nTotalElements;
            h = Th.magnitudes;
            R = abs(V - U);

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
end