classdef ElementFunction < handle
    % ELEMENTFUNCTION Element-level discrete function data.
    %
    %   ElementFunction stores function values and derivatives evaluated at
    %   integration nodes, along with the associated quadrature weights.
    %   This class provides the fundamental data structure for numerical
    %   integration and finite element computations.
    %
    %   The class manages function evaluations at quadrature points and
    %   supports efficient computation of integrals involving functions
    %   and their derivatives. It serves as a bridge between continuous
    %   function representations and discrete numerical computations.
    %
    %   Function derivatives are computed using multi-index notation to
    %   handle arbitrary spatial dimensions and derivative orders
    %   efficiently.
    %
    % Examples:
    %   % Create function data with Gauss quadrature
    %   integrator = approx.integrate.GaussLegendre(5);
    %   funcData = ElementFunction(integrator);
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
    %   weights = funcData.weights;
    %
    % See also:
    %   approx.integrate.Integrator, approx.element.ElementBilinearForm

    properties
        integrator % Integrator object for quadrature
        values % Function values at integration nodes
        derivatives % Function derivatives at integration nodes (cell array)
    end

    properties (Dependent)
        nodes % Integration nodes
        weights % Integration weights
        nPoints % Number of integration points (integer)
        nDerivatives % Number of available derivatives (integer)
    end

    methods
        function obj = ElementFunction(I)
            % ELEMENTFUNCTION Constructor for ElementFunction.
            %
            %   obj = ElementFunction(I) creates a function data object with
            %   the specified integrator. The function values and
            %   derivatives are initially empty and must be set using the
            %   appropriate methods.
            %
            % Inputs:
            %   I - Integrator object for quadrature computations
            %
            % Outputs:
            %   obj - Constructed ElementFunction object

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

            w = obj.integrator.weights;
        end

        function obj = setPoints(obj, varargin)
            % SETPOINTS Set the integration points.
            %
            %   obj = setPoints(obj, varargin) configures the integration
            %   points using the integrator's setPoints method. The
            %   specific arguments depend on the integrator type.
            %
            % Inputs:
            %   obj - The ElementFunction object
            %   varargin - Variable arguments passed to integrator.setPoints
            %
            % Outputs:
            %   obj - The ElementFunction object

            if isempty(varargin)
                return;
            end

            obj.integrator.setPoints(varargin{:});
        end

        function obj = setValues(obj, f, X)
            % SETVALUES Set the function values at specified points.
            %
            %   obj = setValues(obj, f, X) evaluates the function f at the
            %   specified points X and stores the results in the values
            %   property.
            %
            % Inputs:
            %   obj - The ElementFunction object
            %   f - Function object with evaluate method
            %   X - Evaluation nodes (column vector)
            %
            % Outputs:
            %   obj - The ElementFunction object
            
            obj.values = f.evaluate(X);
        end

        function obj = setDerivatives(obj, f, X, k)
            % SETDERIVATIVES Set function derivatives at specified points.
            %
            %   obj = setDerivatives(obj, f, X, k) computes derivatives of
            %   the function @a f up to order @a k at the specified points
            %   @a X using multi-index notation. The derivatives are stored
            %   in the derivatives cell array.
            %
            %   The multi-index approach enables efficient computation of
            %   derivatives in arbitrary spatial dimensions, where each
            %   derivative is characterized by a multi-index specifying
            %   the order of differentiation in each coordinate direction.
            %
            % Inputs:
            %   obj - The ElementFunction object
            %   f - DifferentiableFunction object with derivative method
            %   X - Evaluation nodes (column vector)
            %   k - Maximum derivative order (non-negative integer)
            %
            % Outputs:
            %   obj - The ElementFunction object

            n = size(X, 1);
            
            %< Create multi-index generator for derivative orders
            L = core.linalg.L1MultiIndexer();
            L.setCache(n, n + k);
            
            %< Extract multi-indices (exclude zero-th order)
            M = L.cache(2:end, :) - 1;
            m = size(M, 1);
            
            %< Compute derivatives for each multi-index
            dY = arrayfun(@(i) f.derivative(X, M(i, :)), 1:m, 'Un', 0);
            obj.derivatives = dY;
        end
    end
end