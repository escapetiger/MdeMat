classdef ScaledAffineSpace < approx.space.AffineSpace
    % SCALEDAFFINESPACE Scaled affine transformation of function space.
    %
    %   ScaledAffineSpace extends the affine space with individual scaling
    %   parameters for both the coefficients and the directional component.
    %   Functions in this space have the form:
    %
    %   \f[
    %     f(x) = \varepsilon_u \cdot u(x) + \varepsilon_g \cdot g(x)
    %   \f]
    %
    %   where \f$s_u\f$ are scaling factors applied to the finite element
    %   coefficients, \f$u(x)\f$ is the finite element function, \f$s_g\f$
    %   is the scaling factor for the directional component, and \f$g\f$
    %   is the directional shift.
    %
    % Examples:
    %   % Create scaled affine space with coefficient and direction scaling
    %   coeffScales = [1.5; 2.0; 0.8]; % scale each FE coefficient
    %   directionScale = 0.5;           % scale directional component
    %   scales = [coeffScales; directionScale];
    %   space = ScaledAffineSpace(finiteElement, scales);
    %   
    %   % Evaluate with automatic scaling applied
    %   points = [0, 0.5, 1.0; 0, 0, 0];
    %   coeffs = [1; 1; 1];       % unit coefficients
    %   direction = [1, 1, 1];    % unit direction
    %   values = space.evaluate(points, coeffs, direction);
    %   % Result: values = 1.5*element1 + 2.0*element2 + 0.8*element3 + 0.5*[1,1,1]
    %
    % See also:
    %   approx.space.AffineSpace, approx.space.FunctionSpace

    properties
        scales % Scaling parameters vector [coeffScales; directionScale]
    end

    methods
        function obj = ScaledAffineSpace(element, scales)
            % SCALEDAFFINESPACE Constructor for ScaledAffineSpace.
            %
            %   obj = ScaledAffineSpace(element, scales) creates a scaled
            %   affine space with the specified scaling parameters for
            %   coefficients and directional components.
            %
            % Inputs:
            %   element - Element object for the underlying space
            %   scales - (nDofs+1, 1) scaling parameter vector where:
            %
            % Outputs:
            %   obj - Constructed ScaledAffineSpace object

            if nargin < 2, scales = []; end

            obj@approx.space.AffineSpace(element);
            obj.setScales(scales);
        end

        function obj = setScales(obj, scales)
            % SETSCALES Set the scaling parameters.
            %
            %   obj = setScales(obj, scales) sets the scaling parameters
            %   for the finite element coefficients and directional
            %   component. The scaling vector must have length nDofs + 1.
            %
            % Inputs:
            %   obj - The ScaledAffineSpace object
            %   scales - (nDofs+1, 1) scaling parameter vector
            %
            % Outputs:
            %   obj - The ScaledAffineSpace object

            core.except.assert(isempty(scales) || ...
                length(scales) == (obj.element.nDofs + 1), ...
                'InvalidInput', 'Invalid scaling parameters.');
            obj.scales = scales;
        end

        function Y = evaluate(obj, varargin)
            % EVALUATE Evaluate functions in the scaled affine space.
            %
            %   Y = evaluate(obj, Z) evaluates the scaled affine function
            %   using a combined input matrix, applying scaling to both
            %   coefficients and directional components.
            %
            %   Y = evaluate(obj, Z, G) evaluates with separate directional
            %   component matrix, applying appropriate scaling.
            %
            %   Y = evaluate(obj, X, C, G) evaluates with explicitly
            %   separated components, applying scaling to coefficients
            %   and directional components.
            %
            % Inputs:
            %   obj - The ScaledAffineSpace object
            %   varargin - Input arguments
            %<   Z - Combined matrix [X; C] or [X; C; G]
            %<   X - (nDims, nPoints) evaluation points
            %<   C - (nDofs, nVars) coefficients
            %<   G - (nComponents, nPoints) directional components
            %
            % Outputs:
            %   Y - Scaled function values at evaluation points
            
            sU = obj.scales(1:end-1);
            sG = obj.scales(end);
            E = obj.element;

            switch length(varargin)
                case 1
                    Z = varargin{1};
                    Z = reshape(Z, size(Z, 1), []);
                    n = E.nDims;
                    m = E.nDofs;
                    Z(n+(1:m), :) = sU(:) .* Z(n+(1:m), :);
                    U = E.evaluate(Z(1:(n+m), :));
                    G = sG * Z((n + m + 1):end, :);
                    Y = U + G;
                case 2
                    [Z, G] = varargin{:};
                    n = E.nDims;
                    Z((n+1):end, :) = sU(:) .* Z((n+1):end, :);
                    U = E.evaluate(Z);
                    G = sG * G;
                    Y = U + G;
                case 3
                    [X, C, G] = varargin{:};
                    C = sU(:) .* C;
                    U = E.evaluate(X, C);
                    G = sG * G;
                    Y = U + G;
            end
        end

        function G = direction(obj, Y, varargin)
            % DIRECTION Compute scaled directional component for target
            % values.
            %
            %   G = direction(obj, Y, Z) computes the directional component
            %   needed to reach target values, accounting for scaling
            %   factors applied to both coefficients and the directional
            %   component.
            %
            %   G = direction(obj, Y, X, C) computes the direction using
            %   explicitly separated evaluation points and coefficients,
            %   applying appropriate scaling transformations.
            %
            % Inputs:
            %   obj - The ScaledAffineSpace object
            %   Y - (nComponents, nPoints) target function values
            %   varargin - Input arguments
            %<   Z - (nDims+nDofs, nPoints) combined matrix [X; C]
            %<   X - (nDims, nPoints) evaluation points
            %<   C - (nDofs, nFunctions) finite element coefficients
            %
            % Outputs:
            %   G - (nComponents, nPoints) scaled directional component
            
            sU = obj.scales(1:end-1);
            sG = obj.scales(end);

            switch length(varargin)
                case 1
                    Z = varargin{1};
                    K = obj.element;
                    n = K.nDims;
                    Z((n+1):end, :) = sU(:) .* Z((n+1):end, :);
                    U = obj.element.evaluate(Z);
                    G = (Y - U) / sG;
                case 2
                    [X, C] = varargin{:};
                    C = sU(:) .* C;
                    K = obj.element;
                    U = K.evaluate(X, C);
                    G = (Y - U) / sG;
            end
        end
    end
end