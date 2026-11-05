classdef AffineSpace < approx.space.FunctionSpace
    % AFFINESPACE Affine transformation of function space.
    %
    %   AffineSpace represents an affine transformation of a function
    %   space, allowing for additional degrees of freedom through
    %   directional components. Functions in this space have the form:
    %
    %   \f[
    %     f(x) = u(x) + g(x)
    %   \f]
    %
    %   where \f$u(x)\f$ is a function from the underlying function space
    %   and \f$g\f$ is a directional component.
    %
    % Examples:
    %   % Create affine space from element
    %   space = AffineSpace(finiteElement);
    %   
    %   % Evaluate function with directional component
    %   points = [0, 0.5, 1.0; 0, 0, 0]; % 2D points
    %   coeffs = [1; 0.5; -0.2]; % FE coefficients
    %   direction = [0.1, 0.2, 0.3]; % directional shift
    %   values = space.evaluate(points, coeffs, direction);
    %   
    %   % Compute direction needed to reach target values
    %   targetValues = [1.5, 2.0, 1.8];
    %   requiredDirection = space.direction(targetValues, points, coeffs);
    %
    % See also:
    %   approx.space.FunctionSpace, approx.space.ScaledAffineSpace,
    %   approx.space.MeshSpace
    
    properties (Dependent)
        nDofs   % Number of degrees of freedom
        nPoints % Number of integration points
    end

    methods
        function n = get.nDofs(obj)
            % GET.NDOFS Returns the number of degrees of freedom.

            n = obj.element.nDofs;
        end

        function n = get.nPoints(obj)
            % GET.NPOINTS Returns the number of integration points.

            n = obj.element.volumeData.nPoints;
        end

        function Y = evaluate(obj, varargin)
            % EVALUATE Evaluate functions in the affine space.
            %
            %   Y = evaluate(obj, Z) evaluates the affine function using
            %   a combined input matrix containing points, coefficients,
            %   and directional components.
            %
            %   Y = evaluate(obj, Z, G) evaluates with separate directional
            %   component matrix.
            %
            %   Y = evaluate(obj, X, C, G) evaluates with explicitly
            %   separated points, coefficients, and directional components.
            %
            % Inputs:
            %   obj - The AffineSpace object
            %   varargin - Input arguments
            %<   Z - Combined matrix [X; C] or [X; C; G] where:
            %<       X is (nDims, nPoints) evaluation points
            %<       C is (nDofs, nFunctions) coefficients  
            %<       G is (nComponents, nPoints) directional components
            %<   X - (nDims, nPoints) evaluation points
            %<   C - (nDofs, nFunctions) element coefficients
            %<   G - (nComponents, nPoints) directional components
            %
            % Outputs:
            %   Y - Function values at evaluation points
            
            switch length(varargin)
                case 1
                    Z = varargin{1};
                    Z = reshape(Z, size(Z, 1), []);
                    n = obj.element.nDims;
                    m = obj.element.nDofs;
                    U = obj.element.evaluate(Z(1:(n+m), :));
                    G = Z((n + m + 1):end, :);
                    Y = U + G;
                case 2
                    [Z, G] = varargin{:};
                    U = obj.element.evaluate(Z);
                    Y = U + G;
                case 3
                    [X, C, G] = varargin{:};
                    U = obj.element.evaluate(X, C);
                    Y = U + G;
            end
        end

        function G = direction(obj, Y, varargin)
            % DIRECTION Compute directional component to reach target
            % values.
            %
            %   G = direction(obj, Y, Z) computes the directional component
            %   needed to transform the element approximation to
            %   reach the specified target values.
            %
            %   G = direction(obj, Y, X, C) computes the direction using
            %   explicitly separated evaluation points and coefficients.
            %
            % Inputs:
            %   obj - The AffineSpace object
            %   Y - (nComponents, nPoints) target function values
            %   varargin - Input arguments
            %<   Z - (nDims+nDofs, nPoints) combined matrix [X; C]
            %<   X - (nDims, nPoints) evaluation points
            %<   C - (nDofs, nVars) element coefficients
            %
            % Outputs:
            %   G - (nComponents, nPoints) directional component
            
            U = obj.element.evaluate(varargin{:});
            G = Y - U;
        end
    end
end