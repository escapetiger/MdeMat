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
        derivatives % Function derivatives (nCodims x nPoints x nDerivatives)
        derivativePattern % Multi-index for derivative evaluation 
    end

    properties (Dependent)
        nodes % Integration nodes
        weights % Integration weights
        nPoints % Number of integration points (integer)
        nDerivatives % Number of available derivatives (integer)
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

            obj.setDerivativePattern(order);

            if ~isempty(func)
                obj.setData(func);
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
            % GET.NDERIVATIVES Get the number of available derivatives.

            n = size(obj.derivativePattern, 1);
        end

        function X = get.nodes(obj)
            % GET.NODES Get the integration nodes.

            X = obj.integrator.nodes;
        end

        function w = get.weights(obj)
            % GET.WEIGHTS Get the integration weights.

            w = obj.integrator.weights;
        end

        function obj = setDerivativePattern(obj, order)
            % SETDERIVATIVEPATTERN Set the multi-index pattern for
            % derivative evaluation.
            %
            %   obj = setDerivativePattern(obj, order) generates the
            %   multi-index pattern for derivative evaluation. 
            %
            % Inputs:
            %   obj - The IntegrationData object
            %   order - Derivative order
            %
            % Outputs:
            %   obj - The IntegrationData object

            if order > 0
                nDims = obj.integrat.nDims;
                L = core.linalg.L1MultiIndexer();
                L.setCache(nDims, nDims + order);
                obj.derivativePattern = L.cache(2:end, :) - 1;
            else
                obj.derivativePattern = [];
            end
        end

        function obj = setData(obj, func, varargin)
            % SETDATA Set the integration data.
            %
            %   obj = setData(obj, func, varargin) configures the
            %   integration data by setting integration points, derivative
            %   order, values and derivatives. If varargin = [],
            %   integration points are unchanged.
            %
            % Inputs:
            %   obj - The IntegrationData object
            %   func - Function object
            %   varargin - Arguments passed to integrator.setPoints
            %
            % Outputs:
            %   obj - The IntegrationData object

            if ~isempty(varargin)
                obj.integrator.setPoints(varargin{:});
            end

            obj.setValues(func);
            if ~isempty(obj.derivativePattern)
                obj.setDerivatives(func);
            end
        end

    end

    methods (Access = protected)
        function obj = setValues(obj, func)
            % SETVALUES Set the function values at specified points.
            %
            %   obj = setValues(obj, func) evaluates the function @a func
            %   at the integration points and stores the results in the
            %   values property.
            %
            % Inputs:
            %   obj - The IntegrationData object
            %   func - Function object with evaluate method
            %
            % Outputs:
            %   obj - The IntegrationData object
            
            obj.values = func.evaluate(obj.integrator.nodes);
        end

        function obj = setDerivatives(obj, func)
            % SETDERIVATIVES Set function derivatives at specified points.
            %
            %   obj = setDerivatives(obj, func) computes derivatives of the
            %   function @a func at the integration points using the stored
            %   multi-index pattern.
            %
            %   The multi-index approach enables efficient computation of
            %   derivatives in arbitrary spatial dimensions, where each
            %   derivative is characterized by a multi-index specifying
            %   the order of differentiation in each coordinate direction.
            %
            % Inputs:
            %   obj - The IntegrationData object
            %   func - Function object with derivative method
            %
            % Outputs:
            %   obj - The IntegrationData object

            X = obj.integrator.nodes;
            M = obj.derivativePattern;
            dY = arrayfun(@(i) func.derivative(X, M(i, :)), 1:obj.nDerivatives, 'Un', 0);
            obj.derivatives = cat(3, dY{:});
        end
    end
end