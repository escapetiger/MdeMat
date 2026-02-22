classdef Integrator < handle
    % INTEGRATOR Unified numerical integrator with data storage.
    %
    %   Integrator combines integration point generation, function
    %   evaluation storage, and integration computation in a single
    %   cohesive interface. It manages nodes, weights, function values,
    %   and partial derivatives for efficient numerical integration.
    %
    % See also:
    %   core.geometry.Geometry

    properties
        Rule % Integration rule (approx.integrate.IntegrationRule)
        Nodes % Integration nodes (NDims x NPoints)
        Weights % Integration weights (1 x NPoints)
        Values % Function values (NCodims x NPoints)
        Partials % Partial derivatives (cell of NCodims x NPoints x NPartialsPerOrder)
    end

    properties (Dependent)
        NDims % Number of spatial dimensions
        NPoints % Number of integration points
        NPartials % Number of partial derivatives
        DerivOrder % Maximum derivative order
        NCodims % Number of function components
    end

    methods
        function obj = Integrator(rule)
            % INTEGRATOR Construct an integrator.
            %
            %   obj = Integrator(rule) creates an integrator with
            %   pre-computed integration points.

            arguments
                rule approx.integrate.IntegrationRule
            end

            obj.Rule = rule;
            obj.Nodes = [];
            obj.Weights = [];
            obj.Values = [];
            obj.Partials = [];
        end

        function newObj = copy(obj)
            % COPY Create a deep copy of the integrator.
            %
            %   newObj = copy(obj) creates a new Integrator object with
            %   the same nodes, weights, values, and partial derivatives
            %   as the original object.

            arguments
                obj approx.integrate.Integrator
            end

            newObj = approx.integrate.Integrator(obj.Rule);
            newObj.Nodes = obj.Nodes;
            newObj.Weights = obj.Weights;
            newObj.Values = obj.Values;
            newObj.Partials = obj.Partials;
        end

        function obj = setPoints(obj, geometry, np)
            % SETPOINTS Generate integration points using geometry.
            %
            %   obj = setPoints(obj, geometry, np) sets @a np integration
            %   nodes and weights on the specified geometry @geometry.

            arguments
                obj approx.integrate.Integrator
                geometry core.geometry.Geometry
                np{mustBePositive, mustBeInteger}
            end

            [obj.Nodes, obj.Weights] = obj.Rule.generate(geometry, np);
        end

        function obj = setValues(obj, func, options)
            % SETVALUES Evaluate function at integration points.
            %
            %   obj = setValues(obj, func) evaluates function at the stored
            %   integration nodes and stores the results in the Values
            %   property.
            %
            %   obj = setValues(obj, func, nodes=nodes) evaluates function
            %   at the specified nodes instead of the stored nodes.

            arguments
                obj approx.integrate.Integrator
                func{mustBeFunction}
                options.nodes{mustBeNumeric} = obj.Nodes
            end

            core.except.assert(~isempty(options.nodes), 'InvalidState', ...
                'Integration nodes must be set before evaluating function.');

            if isa(func, 'core.function.Function')
                obj.Values = func.eval(options.nodes);
            else
                obj.Values = func(options.nodes);
            end
        end

        function obj = setPartials(obj, func, order, options)
            % SETPARTIALS Compute partial derivatives at integration
            % points.
            %
            %   obj = setPartials(obj, func, order) computes partial
            %   derivatives up to the specified order at the stored
            %   integration nodes.
            %
            %   obj = setPartials(obj, func, order, nodes=nodes) computes
            %   partial derivatives at the specified nodes instead of the
            %   stored nodes.

            arguments
                obj approx.integrate.Integrator
                func{mustBeFunction}
                order(1, 1) {mustBePositive, mustBeInteger}
                options.nodes{mustBeNumeric} = obj.Nodes
            end

            core.except.assert(~isempty(options.nodes), 'InvalidState', ...
                'Integration nodes must be set before evaluating function.');

            nd = obj.NDims;
            indexer = core.linalg.L1MultiIndexer();
            obj.Partials = cell(1, order);
            for k = 1:order
                indexer.setCache(nd, nd+k, minSum = nd+k);
                M = indexer.Cache - 1;
                dY = cell(1, size(M, 1));
                for j = 1:size(M, 1)
                    dY{j} = func.diff(options.nodes, M(j, :));
                end
                obj.Partials{k} = cat(3, dY{:});
            end
        end

        function obj = addWeightFunction(obj, weightFunc)
            % ADDWEIGHTFUNCTION Apply weight function to integration
            % weights.
            %
            %   obj = addWeightFunction(obj, weightFunc) multiplies the
            %   integration weights by the specified weight function
            %   evaluated at the integration nodes.

            arguments
                obj approx.integrate.Integrator
                weightFunc{mustBeFunction}
            end

            core.except.assert(~isempty(obj.Nodes), 'InvalidState', ...
                'Integration nodes must be set before applying weight function.');

            omega = weightFunc(obj.Nodes);
            omega = reshape(omega, 1, size(obj.Weights, 2));
            obj.Weights = obj.Weights .* omega;
        end

        function I = integrate(obj, func)
            % INTEGRATE Evaluate integral using stored data or function.
            %
            %   I = integrate(obj) integrates using the stored function
            %   values. The Values property must be set first.
            %
            %   I = integrate(obj, func) evaluates the specified function
            %   at the integration nodes and computes the integral.

            arguments
                obj approx.integrate.Integrator
                func{mustBeFunction} = []
            end

            core.except.assert(~isempty(obj.Nodes), 'InvalidState', ...
                'Integration nodes must be set before integration.');

            core.except.assert(~isempty(obj.Weights), 'InvalidState', ...
                'Integration weights must be set before integration.');

            if isempty(func)
                values = obj.Values;
            else
                if isa(func, 'core.function.Function')
                    values = func.eval(obj.Nodes);
                else
                    values = func(obj.Nodes);
                end
            end

            I = values * obj.Weights(:);
        end
    end

    methods
        function n = get.NDims(obj)
            % GET.NDIMS Returns the number of spatial dimensions.

            n = obj.Rule.NDims;
        end

        function n = get.NPoints(obj)
            % GET.NPOINTS Returns the number of integration points.

            if isempty(obj.Weights)
                n = 0;
            else
                n = length(obj.Weights);
            end
        end

        function n = get.NPartials(obj)
            % GET.NPARTIALS Returns the total number of partial derivatives.

            if isempty(obj.Partials)
                n = 0;
            else
                n = sum(cellfun(@(y) size(y, 3), obj.Partials));
            end
        end

        function order = get.DerivOrder(obj)
            % GET.DERIVORDER Returns the maximum derivative order.

            if obj.NPartials == 0 || obj.NDims == 0
                order = 0;
            else
                order = length(obj.Partials);
            end
        end

        function n = get.NCodims(obj)
            % GET.NCODIMS Returns the number of function components.

            if isempty(obj.Values)
                n = 0;
            else
                n = size(obj.Values, 1);
            end
        end
    end
end