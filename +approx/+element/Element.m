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
    %   This abstract base class provides the fundamental interface for
    %   element-based computations, including function evaluation and
    %   access to geometric properties.
    %
    % See also:
    %   core.geometry.Geometry, approx.linear.LinearApproximator

    properties
        geometry % Geometry object
        approximator % LinearApproximator object
    end

    properties (Dependent)
        nDims % Number of spatial dimensions (integer)
        nDofs % Number of degrees of freedom (integer)
    end

    methods
        function obj = Element(geometry, approximator)
            % ELEMENT Constructor for Element.
            %
            %   obj = Element(geometry, approximator) creates an element
            %   with the specified geometry and linear approximator.
            %
            % Inputs:
            %   geometry - Geometry object
            %   approximator - LinearApproximator object
            %
            % Outputs:
            %   obj - Constructed Element object

            obj.geometry = geometry;
            obj.approximator = approximator;
        end

        function n = get.nDims(obj)
            % GET.NDIMS Get the number of spatial dimensions.

            n = obj.approximator.basis.nDims;
        end

        function n = get.nDofs(obj)
            % GET.NDOFS Get the number of degrees of freedom.

            n = obj.approximator.nDofs;
        end

        function Y = evaluate(obj, varargin)
            % EVALUATE Evaluate a function over the element.
            %
            %   Y = evaluate(obj, Z) evaluates the basis functions at
            %   points specified in Z = [X; C], where X contains evaluation
            %   points and C contains coefficients.
            %
            %   Y = evaluate(obj, X, C) evaluates the finite element
            %   function defined by coefficients C at the evaluation
            %   points X.
            %
            % Inputs:
            %   obj - The Element object
            %   varargin - Variable input arguments:
            %<    Z - Combined matrix [X; C] (for single argument case)
            %<    X - Evaluation points (nDims x nPts matrix)
            %<    C - Coefficient vectors (nDofs x nVars matrix)
            %
            % Outputs:
            %   Y - Function values at evaluation points:
            %<     * If nargin = 2: (1 x nPts) vector
            %<     * If nargin = 3: (nVars x nPts) matrix

            Y = obj.approximator.evaluate(varargin{:});
        end
    end
end