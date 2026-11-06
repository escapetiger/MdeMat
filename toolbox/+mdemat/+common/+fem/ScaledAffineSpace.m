classdef ScaledAffineSpace < mdemat.common.fem.AffineSpace
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
    % See also:
    %   mdemat.common.fem.AffineSpace
    
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
            %   obj - Constructed ScaledAffineSpace object
            
            arguments
                element mdemat.common.fem.Element
                scales (:,1) {mustBeNumeric, mustBeFinite}
            end
            
            obj@mdemat.common.fem.AffineSpace(element);
            obj.setScales(scales);
        end
        
        function obj = setScales(obj, scales)
            % SETSCALES Set the scaling parameters.
            %
            %   obj = setScales(obj, scales) sets the scaling parameters
            %   for the finite element coefficients and directional
            %   component. The scaling vector must have length nDofs + 1.
            
            arguments
                obj mdemat.common.fem.ScaledAffineSpace
                scales (:,1) {mustBeNumeric, mustBeFinite}
            end
            
            assert(all(isempty(scales) || ...
                length(scales) == (obj.Element.NDofs + 1), 'all'), ...
                'mdemat:common:fem:ScaledAffineSpace:InvalidInput', '%s', 'Invalid scaling parameters.');
            if isempty(scales)
                scales = ones(obj.Element.NDofs + 1, 1);
            end
            obj.scales = scales;
        end
        
        function Y = eval(obj, X, C, G)
            % EVAL Evaluate functions in the scaled affine space.
            %
            %   Y = eval(obj, X, C, G) evaluates with explicitly
            %   separated components, applying scaling to coefficients
            %   and directional components.
            
            arguments
                obj mdemat.common.fem.ScaledAffineSpace
                X {mustBeNumeric}
                C {mustBeNumeric}
                G {mustBeNumeric}
            end
            
            sU = obj.scales(1:end-1);
            sG = obj.scales(end);
            C = sU(:) .* C;
            U = obj.Element.eval(X, C);
            G = sG * G;
            Y = U + G;
        end
        
        function G = direction(obj, Y, X, C)
            % DIRECTION Compute scaled directional component for target
            % values.
            %
            %   G = direction(obj, Y, X, C) computes the direction using
            %   explicitly separated evaluation points and coefficients,
            %   applying appropriate scaling transformations.

            arguments
                obj mdemat.common.fem.ScaledAffineSpace
                Y {mustBeNumeric}
                X {mustBeNumeric}
                C {mustBeNumeric}
            end
            
            sU = obj.scales(1:end-1);
            sG = obj.scales(end);
            C = sU(:) .* C;
            U = obj.Element.eval(X, C);
            assert(all(sG ~= 0, 'all'), 'mdemat:common:fem:ScaledAffineSpace:ZeroDirectionScale', ...
                '%s', 'The directional scale must be nonzero to recover a direction.');
            G = (Y - U) / sG;
        end
    end
end