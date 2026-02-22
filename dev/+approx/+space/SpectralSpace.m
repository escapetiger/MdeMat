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
    % See also:
    %   approx.element.Element
    
    properties
        Element % Element object defining the function space
    end
    
    properties (Dependent)
        NDofs % Number of degrees of freedom (integer)
        NDims % Number of spatial dimensions (integer)
    end
    
    methods
        function obj = SpectralSpace(element)
            % SPECTRALSPACE Constructor for SpectralSpace.
            %
            %   obj = SpectralSpace(element) creates a spectral space
            %   using the specified element as the basis for the
            %   finite-dimensional function space.

            arguments
                element approx.element.Element
            end
            
            obj.Element = element;
        end
        
        function n = get.NDofs(obj)
            % GET.NDOFS Get the number of degrees of freedom.
            
            n = obj.Element.NDofs;
        end
        
        function n = get.NDims(obj)
            % GET.NDIMS Get the number of spatial dimensions.
            
            n = obj.Element.NDims;
        end
        
        function Y = eval(obj, X, C)
            % EVAL Evaluate functions in the spectral space.
            %
            %   Y = eval(obj, X, C) evaluates functions with
            %   coefficient matrix @a C at evaluation points @a X.

            arguments
                obj approx.space.SpectralSpace
                X {mustBeNumeric}
                C {mustBeNumeric}
            end
            
            Y = obj.Element.eval(X, C);
        end

        function U = project(obj, C)
            % PROJECT Project coefficients onto the spectral space.
            %
            %   U = project(obj, C) projects the coefficient matrix @a C
            %   onto the spectral space defined by the element.

            arguments
                obj approx.space.SpectralSpace
                C {mustBeNumeric}
            end

            U = obj.Element.project(C);
        end
    end
end