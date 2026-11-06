classdef Integrator < handle
    % INTEGRATOR Unified numerical integrator with data storage.
    %
    %   Integrator combines integration point generation, function
    %   evaluation storage, and integration computation in a single
    %   cohesive interface. It manages nodes, weights, function values,
    %   and partial derivatives for efficient numerical integration.
    %
    % See also:
    %   mdemat.data.geometry.Geometry

    properties
        Rule % Integration rule (mdemat.common.quadrature.IntegrationRule)
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
                rule mdemat.common.quadrature.IntegrationRule
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
                obj mdemat.common.quadrature.Integrator
            end

            newObj = mdemat.common.quadrature.Integrator(obj.Rule);
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
                obj mdemat.common.quadrature.Integrator
                geometry mdemat.data.geometry.Geometry
                np{mustBePositive, mustBeInteger}
            end

            [obj.Nodes, obj.Weights] = obj.Rule.generate(geometry, np);
            obj.Values = [];
            obj.Partials = [];
        end

        function obj = setBoundaryPoints(obj, geometry, np)
            % SETBOUNDARYPOINTS Set physical face quadrature, including traces.
            %   A collapsed orthotope describes a face with its intrinsic
            %   surface measure; setPoints instead uses ambient volume.
            arguments
                obj mdemat.common.quadrature.Integrator
                geometry mdemat.data.geometry.Orthotope
                np {mustBeVector, mustBePositive, mustBeInteger}
            end
            assert(isa(obj.Rule, 'mdemat.common.quadrature.OrthotopeRule'), ...
                'MdeMat:Quadrature:UnsupportedBoundaryRule', ...
                'Axis-aligned face quadrature requires an orthotope rule.');
            [obj.Nodes, obj.Weights] = obj.Rule.generateOnBoundary(geometry, np);
            obj.Values = [];
            obj.Partials = [];
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
                obj mdemat.common.quadrature.Integrator
                func{mustBeFunction}
                options.nodes{mustBeNumeric} = obj.Nodes
            end

            assert(all(~isempty(options.nodes), 'all'), 'mdemat:common:quadrature:Integrator:InvalidState', ...
                '%s', 'Integration nodes must be set before evaluating function.');

            if isa(func, 'mdemat.data.function.Function')
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
                obj mdemat.common.quadrature.Integrator
                func{mustBeFunction}
                order(1, 1) {mustBePositive, mustBeInteger}
                options.nodes{mustBeNumeric} = obj.Nodes
            end

            assert(all(~isempty(options.nodes), 'all'), 'mdemat:common:quadrature:Integrator:InvalidState', ...
                '%s', 'Integration nodes must be set before evaluating function.');

            nd = obj.NDims;
            indexer = mdemat.data.index.L1MultiIndexer();
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
                obj mdemat.common.quadrature.Integrator
                weightFunc{mustBeFunction}
            end

            assert(all(~isempty(obj.Nodes), 'all'), 'mdemat:common:quadrature:Integrator:InvalidState', ...
                '%s', 'Integration nodes must be set before applying weight function.');

            if isa(weightFunc, 'mdemat.data.function.Function')
                omega = weightFunc.eval(obj.Nodes);
            else
                omega = weightFunc(obj.Nodes);
            end
            if isscalar(omega)
                omega = repmat(omega, 1, obj.NPoints);
            end
            assert(all(isnumeric(omega) && numel(omega) == obj.NPoints, 'all'), ...
                'mdemat:common:quadrature:Integrator:InvalidWeightShape', '%s', 'Weight values must be scalar or one per node.');
            omega = reshape(omega, 1, obj.NPoints);
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
                obj mdemat.common.quadrature.Integrator
                func = []
            end

            assert(all(~isempty(obj.Nodes), 'all'), 'mdemat:common:quadrature:Integrator:InvalidState', ...
                '%s', 'Integration nodes must be set before integration.');

            assert(all(~isempty(obj.Weights), 'all'), 'mdemat:common:quadrature:Integrator:InvalidState', ...
                '%s', 'Integration weights must be set before integration.');

            if ~isempty(func)
                mustBeFunction(func);
            end
            if isempty(func)
                values = obj.Values;
            else
                if isa(func, 'mdemat.data.function.Function')
                    values = func.eval(obj.Nodes);
                else
                    values = func(obj.Nodes);
                end
            end

            assert(all(~isempty(values), 'all'), 'mdemat:common:quadrature:Integrator:InvalidState', ...
                '%s', 'No function values have been supplied.');
            if isscalar(values)
                values = repmat(values, 1, obj.NPoints);
            end
            assert(all(size(values, 2) == obj.NPoints, 'all'), ...
                'mdemat:common:quadrature:Integrator:InvalidValueShape', '%s', 'Function values must have one column per node.');
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
