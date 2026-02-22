classdef Element < handle
    % ELEMENT Base class for element representations.
    %
    %   Element encapsulates geometry information and function space used
    %   in Galerkin methods. The main idea of Galerkin methods is to
    %   approximate an infinite-dimensional function space by a
    %   finite-dimensional function space, called the "ansatz space".
    %
    %   In the spectral method (SPM), the ansatz space consists of
    %   high-order polynomials on an element. In the finite element method
    %   (FEM), the ansatz space comprises piecewise low-order polynomials
    %   over a mesh of elements. Element serves as a building block for
    %   both SPM and FEM to define the ansatz space.
    %
    % See also:
    %   core.geometry.Geometry, approx.linear.LinearApproximator

    properties
        Geometry % Geometry object
        Approximator % LinearApproximator object
    end

    properties (Dependent)
        NDims % Number of spatial dimensions (integer)
        NDofs % Number of degrees of freedom (integer)
    end

    methods
        function obj = Element(geometry, approximator)
            % ELEMENT Constructor for Element.
            %
            %   obj = Element(geometry, approximator) creates an element
            %   with the specified @a geometry and @a approximator.

            arguments
                geometry core.geometry.Geometry
                approximator approx.linear.LinearApproximator
            end

            obj.Geometry = geometry;
            obj.Approximator = approximator;
        end

        function n = get.NDims(obj)
            % GET.NDIMS Get the number of spatial dimensions.

            n = obj.Approximator.Basis.NDims;
        end

        function n = get.NDofs(obj)
            % GET.NDOFS Get the number of degrees of freedom.

            n = obj.Approximator.NDofs;
        end

        function Y = eval(obj, X, C)
            % EVAL Evaluate a function over the element.
            %
            %   Y = eval(obj, X, C) evaluates the finite element
            %   function defined by coefficients @a C at the evaluation
            %   points @a X.

            arguments
                obj approx.element.Element
                X {mustBeNumeric, mustBeNonempty}
                C {mustBeNumeric, mustBeNonempty}
            end

            Y = obj.Approximator.eval(X, C);
        end
        
        function Y = diff(obj, X, C, order)
            % DIFF Differentiate a function over the element.
            %
            %   Y = diff(obj, X, C, order) differentiates the
            %   finite element function defined by coefficients @a C at the
            %   evaluation points @a X with derivative @a order.
            
            arguments
                obj approx.element.Element
                X {mustBeNumeric, mustBeNonempty}
                C {mustBeNumeric, mustBeNonempty}
                order (1, :) {mustBeNumeric, mustBeNonempty}
            end
            
            Y = obj.Approximator.diff(X, C, order);
        end

        function Y = fit(obj, F)
            % FIT Fit data over the element.
            %
            %   U = fit(obj, F) fits data over the element by
            %   solving the fitting system \f$D*U = F\f$, where \f$D\f$ is
            %   the design matrix. Handles well-posed and overdetermined
            %   systems using least squares when necessary.
            
            arguments
                obj approx.element.Element
                F {mustBeNumeric, mustBeNonempty}
            end

            Y = obj.Approximator.fit(F);
        end

        function Y = project(obj, F)
            % PROJECT Project data over the element.
            %
            %   U = project(obj, F) projects data over the element by
            %   solving the fitting system \f$D*U = F\f$, where \f$D\f$ is
            %   the design matrix. Handles well-posed and overdetermined
            %   systems using least squares when necessary.
            
            arguments
                obj approx.element.Element
                F {mustBeNumeric, mustBeNonempty}
            end

            Y = obj.Approximator.project(F);
        end
    end
end