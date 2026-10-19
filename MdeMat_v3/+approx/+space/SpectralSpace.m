classdef SpectralSpace < handle
    % SPECTRALSPACE Function space based on a single element.
    %
    %   SpectralSpace represents a finite-dimensional function space
    %   defined by a single element. It serves as a building block for
    %   more complex function space constructions through direct sum
    %   operations.
    %
    %   Functions in this space have the form:
    %   \f[
    %       u_h(x) = \sum_{i=1}^n c_i \phi_i(x)
    %   \f]
    %   where \f$\phi_i\f$ are basis functions from the element and
    %   \f$c_i\f$ are coefficients.
    %
    % Examples:
    %   % Create spectral space from element
    %   element = approx.element.C0Element(geometry, approximator, integrator);
    %   space = approx.space.SpectralSpace(element);
    %   
    %   % Evaluate function at points
    %   points = [0, 0.5, 1.0; 0, 0, 0];  % 2D points
    %   coeffs = [1; 0.5; -0.2];          % Coefficients
    %   values = space.evaluate(points, coeffs);
    %
    % See also:
    %   approx.space.DirectSum, approx.element.Element

    properties
        element % Element object defining the function space
    end

    properties (Dependent)
        nDofs % Number of degrees of freedom (integer)
        nDims % Number of spatial dimensions (integer)
    end

    methods
        function obj = SpectralSpace(element)
            % SPECTRALSPACE Constructor for SpectralSpace.
            %
            %   obj = SpectralSpace(element) creates a spectral space
            %   using the specified element as the basis for the
            %   finite-dimensional function space.
            %
            % Inputs:
            %   element - Element object
            %
            % Outputs:
            %   obj - Constructed SpectralSpace object
            
            obj.element = element;
        end

        function n = get.nDofs(obj)
            % GET.NDOFS Get the number of degrees of freedom.

            n = obj.element.nDofs;
        end

        function n = get.nDims(obj)
            % GET.NDIMS Get the number of spatial dimensions.

            n = obj.element.nDims;
        end

        function Y = evaluate(obj, X, C)
            % EVALUATE Evaluate functions in the spectral space.
            %
            %   Y = evaluate(obj, X, C) evaluates functions with
            %   coefficient matrix @a C at evaluation points @a X.
            %
            % Inputs:
            %   obj - The SpectralSpace object
            %   X - Evaluation points (nDims x nPoints matrix)
            %   C - Coefficient matrix (nDofs x nComponents matrix)
            %
            % Outputs:
            %   Y - Function values at evaluation points (nComponents x nPoints matrix)
            
            Y = obj.element.evaluate(X, C);
        end
    end
end