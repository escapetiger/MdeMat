classdef IntegrationData < handle
    % INTEGRATIONDATA Data collection for integration.
    %
    %   IntegrationData stores function values and derivatives evaluated at
    %   integration nodes, along with the associated integration weights.
    %
    %   The class manages function evaluations at integration points and
    %   supports efficient computation of integrals involving functions
    %   and their derivatives. It serves as a bridge between continuous
    %   function representations and discrete numerical computations.
    %
    %   Function derivatives are computed using multi-index notation to
    %   handle arbitrary spatial dimensions and derivative orders
    %   efficiently.
    %
    % See also:
    %   approx.integrate.Integrator

    properties
        integrator % Integrator object
        values % Function values (nCodims x nPoints)
        derivatives % Function derivatives (cell of nCodims x nPoints x nDerivatives)
    end

    properties (Dependent)
        nodes % Integration nodes
        weights % Integration weights
        nPoints % Number of integration points (integer)
        nDerivatives % Number of all derivatives (integer)
        derivativeOrder % Derivative order
    end

    methods
        function obj = IntegrationData(integrator, order, func)
            % INTEGRATIONDATA Constructor for IntegrationData.
            %
            %   obj = IntegrationData(integrator, order) creates an
            %   integration data object with the specified integrator and
            %   derivative order. Values and 
            %
            %   obj = IntegrationData(integrator, order) creates an
            %   integration data object with the specified integrator and
            %   derivative order.
            %
            % Inputs:
            %   integrator - Integrator object
            %   func - Function object
            %   order - Derivative order
            %
            % Outputs:
            %   obj - Constructed IntegrationData object

            if nargin < 3
                func = [];
            end

            obj.integrator = integrator;
            if ~isempty(func)
                obj.setValues(func);
                obj.setDerivatives(func, order);
            else
                obj.values = [];
                obj.derivatives = {};
            end
        end

        function n = get.nPoints(obj)
            % GET.NPOINTS Get the number of integration points.

            n = length(obj.integrator.weights);
        end

        function n = get.nDerivatives(obj)
            % GET.NDERIVATIVES Get the number of all derivatives.

            n = sum(cellfun(@(y) size(y, 3), obj.derivatives));
        end

        function n = get.derivativeOrder(obj)
            % GET.DERIVATIVEORDER Get the derivative order.

            n = length(obj.derivatives);
        end

        function X = get.nodes(obj)
            % GET.NODES Get the integration nodes.

            X = obj.integrator.nodes;
        end

        function w = get.weights(obj)
            % GET.WEIGHTS Get the integration weights.

            w = obj.integrator.weights;
        end

        function obj = setValues(obj, func, nodes)
            % SETVALUES Set the function values at specified points.
            %
            %   obj = setValues(obj, func) evaluates the function @a func
            %   at the integration points and stores the results in the
            %   values property.
            %
            %   obj = setValues(obj, func, nodes) evaluates the function @a
            %   func at the specified integration points @a nodes and
            %   stores the results in the values property.
            %
            % Inputs:
            %   obj - The IntegrationData object
            %   func - Function object with evaluate method
            %   nodes - External integration nodes
            %
            % Outputs:
            %   obj - The IntegrationData object
            
            if nargin < 3 || isempty(nodes)
                nodes = obj.integrator.nodes;
            end
            
            obj.values = func.evaluate(nodes);
        end

        function obj = setDerivatives(obj, func, order, nodes)
            % SETDERIVATIVES Set function derivatives at specified points.
            %
            %   obj = setDerivatives(obj, func, order) computes derivatives
            %   of the function @a func up to order @a order at the
            %   integration points.
            %
            %   obj = setDerivatives(obj, func, order, nodes) computes
            %   derivatives of the function @a func up to order @a order at
            %   the specified integration points @a nodes.
            %
            % Inputs:
            %   obj - The IntegrationData object
            %   func - Function object with differentiate method
            %   order - Derivative order
            %   nodes - External integration nodes
            %
            % Outputs:
            %   obj - The IntegrationData object

            if nargin < 4  || isempty(nodes)
                nodes = obj.integrator.nodes;
            end

            d = obj.integrator.nDims;
            X = nodes;
            L = core.linalg.L1MultiIndexer();
            for k = 1:order
                L.setCache(d, d + k, d + k);
                M = L.cache - 1;
                dY = arrayfun(@(i) func.differentiate(X, M(i, :)), 1:size(M, 1), 'Un', 0);
                obj.derivatives{k} = cat(3, dY{:});
            end
        end
    end
end