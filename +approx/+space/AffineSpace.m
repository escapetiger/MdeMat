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
    % See also:
    %   approx.space.ScaledAffineSpace
    
    properties (Dependent)
        NDofs   % Number of degrees of freedom
        NPoints % Number of integration points
        NDims % Number of dimensions
    end
    
    properties
        Element % Element object
    end
    
    methods
        function obj = AffineSpace(element)
            % AFFINESPACE Constructor for AffineSpace.
            %
            %   obj = AffineSpace(element) creates an affine space using
            %   the specified element representation.

            arguments
                element approx.element.Element
            end
            
            obj.Element = element;
        end
        
        function n = get.NDofs(obj)
            % GET.NDOFS Returns the number of degrees of freedom.
            
            n = obj.Element.NDofs;
        end
        
        function n = get.NPoints(obj)
            % GET.NPOINTS Returns the number of integration points.
            
            n = obj.Element.Volume.NPoints;
        end
        
        function Y = eval(obj, X, C, G)
            % EVAL Evaluate functions in the affine space.
            %
            %   Y = eval(obj, X, C, G) evaluates with explicitly
            %   separated points, coefficients, and directional components.
            
            arguments
                obj approx.space.AffineSpace
                X {mustBeNumeric}
                C {mustBeNumeric}
                G {mustBeNumeric}
            end
            
            U = obj.Element.eval(X, C);
            Y = U + G;
        end
        
        function G = direction(obj, Y, X, C)
            % DIRECTION Compute directional component to reach target
            % values.
            %
            %   G = direction(obj, Y, X, C) computes the direction using
            %   explicitly separated evaluation points and coefficients.
            
            arguments
                obj approx.space.AffineSpace
                Y {mustBeNumeric}
                X {mustBeNumeric}
                C {mustBeNumeric}
            end
            
            U = obj.Element.eval(X, C);
            G = Y - U;
        end
    end
end