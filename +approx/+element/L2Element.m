classdef L2Element < approx.element.Element
    % L2ELEMENT \f$L^2\f$ element with volume and flux integration.
    %
    %   L2Element implements a discontinuous finite element based on the
    %   \f$L^2\f$ function space. This class extends the base Element class
    %   by adding volume and flux integration capabilities through
    %   precomputed basis values and derivatives at integration points.
    % 
    %   L2Element assumes discontinuity occurs only at element boundaries,
    %   making them suitable for discontinuous Galerkin (DG) methods and
    %   finite volume schemes. The element supports both modal (orthogonal)
    %   and nodal (interpolation) basis functions.
    %
    % See also:
    %   approx.integrate.IntegrationData, approx.element.Element

    properties (Constant)
        TYPE = 'l2' % Element type identifier
    end

    properties
        volume % Volume integration data (IntegrationData)
        flux % Flux integration data (IntegrationData array)
    end

    properties (Dependent)
        nFluxes % Number of fluxes (integer)
    end

    methods
        function obj = L2Element(G, A, VI, FI)
            % L2Element Constructor for L2Element.
            %
            %   obj = L2Element(G, A, VI, FI) creates a L2 element with
            %   both volume and flux integration capabilities using the
            %   specified components.
            %
            % Inputs:
            %   G - Geometry object
            %   A - LinearApproximator object 
            %   VI - Integrator object
            %   FI - Array of Integrator objects
            %
            % Outputs:
            %   obj - Constructed L2Element object

            obj@approx.element.Element(G, A);
            obj.volume = approx.integrate.IntegrationData(VI, 1, A.basis);
            obj.flux = arrayfun(@(I) approx.integrate.IntegrationData(I, 0, A.basis), FI);
        end

        function n = get.nFluxes(obj)
            % GET.NFLUXES Get the number of fluxes.

            n = length(obj.flux);
        end

        function Y = differentiate(obj, varargin)
            % DIFFERENTIATE Differentiate a function over the element.
            %
            %   Y = differentiate(obj, Z, order) differentiates the basis
            %   functions at points specified in Z = [X; C], where X
            %   contains evaluation points and C contains coefficients.
            %
            %   Y = differentiate(obj, X, C, order) differentiates the
            %   finite element function defined by coefficients C at the
            %   evaluation points X.
            %
            % Inputs:
            %   obj - The Element object
            %   varargin - Variable input arguments:
            %<    Z - Combined matrix [X; C] (for single argument case)
            %<    X - Evaluation points (nDims x nPts matrix)
            %<    C - Coefficient vectors (nDofs x nVars matrix)
            %<    order - Derivative order (nDims x 1 vector)
            %
            % Outputs:
            %   Y - Function values at evaluation points:
            %<     * If nargin = 2: (1 x nPts) vector
            %<     * If nargin = 3: (nVars x nPts) matrix

            Y = obj.approximator.differentiate(varargin{:});
        end
    end
end