classdef FiniteElement < handle
    % FINITEELEMENT Abstract base class for finite element representations.
    %
    %   FiniteElement provides the foundation for finite element methods
    %   with domain geometry, function space, and integration capabilities.
    %   This abstract base class defines the common interface for all
    %   finite element implementations in the finite element method
    %   framework.
    %
    %   The class encapsulates the essential components of a finite
    %   element: the geometric domain and the function space projector. It
    %   provides access to geometric properties and function evaluation
    %   capabilities that are inherited by all concrete finite element
    %   types.
    %
    %   Finite elements discretize continuous function spaces by defining
    %   basis functions over geometric domains and providing numerical
    %   integration capabilities for computing integrals in weak
    %   formulations.
    %
    % See also:
    %   fem.element.C0FiniteElement, fem.element.L2FiniteElement,
    %   core.geometry.Geometry, approx.project.Projector

    properties
        geometry % Geometry object defining the element domain
        projector % Function space projector for basis functions
    end
    
    properties (Dependent)
        nDims % Number of spatial dimensions
        nDofs % Number of degrees of freedom
        outerNormals % Outer normal vectors at geometry boundaries
    end

    methods
        function obj = FiniteElement(K, P)
            % FINITEELEMENT Constructor for FiniteElement.
            %
            %   obj = FiniteElement(K, P) creates a finite element with the
            %   specified geometry and function space projector. This
            %   constructor is called by concrete subclasses to initialize
            %   the base class functionality.
            %
            % Inputs:
            %   K - Element geometry object (Geometry)
            %   P - Function space projector (Projector)
            %
            % Outputs:
            %   obj - Constructed FiniteElement object
            
            obj.geometry = K;
            obj.projector = P;
        end

        function n = get.nDims(obj)
            % GET.NDIMS Get the spatial dimension of the function domain.
            %
            %   n = get.nDims(obj) returns the number of spatial dimensions
            %   of the function domain from the underlying basis functions.
            %
            % Inputs:
            %   obj - The FiniteElement object
            %
            % Outputs:
            %   n - Number of spatial dimensions (positive integer)
            
            n = obj.projector.basis.nDims;
        end

        function n = get.nDofs(obj)
            % GET.NDOFS Get the number of degrees of freedom.
            %
            %   n = get.nDofs(obj) returns the dimension of the finite
            %   element function space, which equals the number of basis
            %   functions or degrees of freedom.
            %
            % Inputs:
            %   obj - The FiniteElement object
            %
            % Outputs:
            %   n - Number of degrees of freedom (positive integer)
            
            n = obj.projector.nDofs;
        end

        function N = get.outerNormals(obj)
            % GET.OUTERNORMALS Get outer normal vectors at geometry
            % boundaries.
            %
            %   N = get.outerNormals(obj) returns the outward-pointing unit
            %   normal vectors at the boundaries of the element geometry.
            %   These are essential for flux computations and boundary
            %   condition enforcement.
            %
            % Inputs:
            %   obj - The FiniteElement object
            %
            % Outputs:
            %   N - Outer normal vectors at boundary points (matrix)
            
            N = obj.geometry.outerNormals();
        end

        function Y = evaluate(obj, varargin)
            % EVALUATE Evaluate functions in the finite element space.
            %
            %   Y = evaluate(obj, Z) evaluates the basis functions at
            %   points specified in the combined matrix Z = [X; C], where X
            %   contains evaluation points and C contains coefficients.
            %
            %   Y = evaluate(obj, X, C) evaluates the finite element
            %   function defined by coefficients C at the evaluation points
            %   X.
            %
            % Inputs:
            %   obj - The FiniteElement object
            %   varargin - Input arguments
            %<   Z - Combined matrix [X; C] where X are points and C are coefficients
            %<   X - Evaluation points (nDims x nPoints matrix)
            %<   C - Coefficient vectors (nDofs x nFunctions matrix)
            %
            % Outputs:
            %   Y - Function values at evaluation points
            %<       * If nargin = 2: (1 x nPoints) vector of basis evaluations
            %<       * If nargin = 3: (nFunctions x nPoints) matrix of function values
            
            Y = obj.projector.evaluate(varargin{:});
        end
    end
end