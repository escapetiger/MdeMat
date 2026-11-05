classdef FunctionData < handle
    % FUNCTIONDATA A collection of discrete data for a function.
    %
    %   FunctionData stores function values and derivatives evaluated at
    %   integration nodes, along with the associated quadrature weights.
    %   This class provides the fundamental data structure for numerical
    %   integration and finite element computations.
    %
    %   The class manages function evaluations at quadrature points and
    %   supports efficient computation of integrals involving functions
    %   and their derivatives. It serves as a bridge between continuous
    %   function representations and discrete numerical computations.
    %
    % Examples:
    %   % Create function data with Gauss quadrature
    %   integrator = approx.integrate.GaussLegendre(5);
    %   funcData = FunctionData(integrator);
    %
    %   % Set integration points and evaluate function
    %   funcData.setPoints(-1, 1);
    %   f = @(x) sin(pi*x);
    %   funcData.setValues(f, funcData.nodes);
    %
    %   % Compute derivatives
    %   df = DifferentiableFunction(f);
    %   funcData.setDerivatives(df, funcData.nodes, 2);
    %
    %   % Access computed data
    %   values = funcData.values;
    %   firstDeriv = funcData.derivatives{1};
    %
    % See also:
    %   approx.integrate.Integrator, approx.fem.data.BilinearFormData

    properties
        integrator % Integrator object for quadrature
        values = [] % Function values at integration nodes
        derivatives = {} % Function derivatives at integration nodes (cell array)
    end

    properties (Dependent)
        nodes % Integration nodes (dependent property)
        weights % Integration weights (dependent property)
        nPoints % Number of integration points (dependent property)
        nDerivatives % Number of available derivatives (dependent property)
    end

    methods
        function obj = FunctionData(I)
            % FUNCTIONDATA Constructor for FunctionData.
            %
            %   obj = FunctionData(I) creates a function data object with
            %   the specified integrator. The function values and derivatives
            %   are initially empty and must be set using the appropriate
            %   methods.
            %
            % Inputs:
            %   I - Integrator object for quadrature computations
            %
            % Outputs:
            %   obj - Constructed FunctionData object

            obj.integrator = I;
            obj.values = [];
            obj.derivatives = {};
        end

        function n = get.nPoints(obj)
            % GET.NPOINTS Get the number of integration points.

            n = length(obj.integrator.weights);
        end

        function n = get.nDerivatives(obj)
            % GET.NDERIVATIVES Get the number of available derivatives.

            n = length(obj.derivatives);
        end

        function X = get.nodes(obj)
            % GET.NODES Get the integration nodes.

            X = obj.integrator.nodes;
        end

        function w = get.weights(obj)
            % GET.WEIGHTS Get the integration weights.
            %
            %   w = get.weights(obj) returns the integration weights from
            %   the associated integrator object.
            %
            % Inputs:
            %   obj - The FunctionData object
            %
            % Outputs:
            %   w - Integration weights (column vector)

            w = obj.integrator.weights;
        end

        function setPoints(obj, varargin)
            % SETPOINTS Set the integration points.
            %
            %   setPoints(obj, varargin) configures the integration points
            %   using the integrator's setPoints method. The specific
            %   arguments depend on the integrator type.
            %
            % Inputs:
            %   obj - The FunctionData object
            %   varargin - Variable arguments passed to integrator.setPoints
            %
            % Outputs:
            %   NULL

            obj.integrator.setPoints(varargin{:});
        end

        function setValues(obj, f, X)
            % SETVALUES Set the function values at specified points.
            %
            %   setValues(obj, f, X) evaluates the function f at the
            %   specified points X and stores the results in the values
            %   property.
            %
            % Inputs:
            %   obj - The FunctionData object
            %   f - Function object with evaluate method
            %   X - Evaluation nodes (column vector)
            %
            % Outputs:
            %   NULL
            
            Y = f.evaluate(X);
            obj.values = Y;
        end

        function setDerivatives(obj, f, X, k)
            % SETDERIVATIVES Set the function derivatives at specified
            % points.
            %
            %   setDerivatives(obj, f, X, k) computes derivatives of the
            %   function f up to order k at the specified points X using
            %   multi-index notation. The derivatives are stored in the
            %   derivatives cell array.
            %
            % Inputs:
            %   obj - The FunctionData object
            %   f - DifferentiableFunction object with derivative method
            %   X - Evaluation nodes (column vector)
            %   k - Maximum derivative order (non-negative integer)
            %
            % Outputs:
            %   NULL

            n = size(X, 1);
            
            %< Create multi-index generator for derivative orders
            L = core.linalg.L1MultiIndexer();
            L.setCache(n, n + k);
            
            %< Extract multi-indices (exclude zero-th order)
            M = L.cache(2:end, :) - 1;
            m = size(M, 1);
            
            %< Compute derivatives for each multi-index
            dY = arrayfun(@(i) f.derivative(X, M(i, :)), 1:m, 'UniformOutput', false);
            obj.derivatives = dY;
        end
    end
end