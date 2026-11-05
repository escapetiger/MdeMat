classdef MeshSpace < fem.space.FiniteElementSpace
    % MESHSPACE Finite element function space defined on a geometric mesh.
    %
    %   MeshSpace extends finite element spaces to operate on geometric
    %   meshes, providing capabilities for function evaluation, projection,
    %   and integration across multiple elements in a discretized domain.
    %
    %   This class handles the mapping between reference elements and
    %   physical elements in the mesh, allowing for efficient computation
    %   of global operations such as L2 projections and weighted averages.
    %   It serves as the foundation for finite element computations on
    %   complex geometries.
    %
    % Examples:
    %   % Create mesh space from finite element and mesh
    %   space = MeshSpace(finiteElement, mesh);
    %   
    %   % Project analytical function onto the mesh space
    %   analyticFunc = @(x) sin(pi*x(1,:)) .* cos(pi*x(2,:));
    %   coefficients = space.project(analyticFunc);
    %   
    %   % Evaluate projected function at specific points
    %   evalPoints = [0.25, 0.75; 0.5, 0.5];
    %   values = space.evaluate(evalPoints, coefficients);
    %   
    %   % Compute weighted average of a function over the domain
    %   quadraticFunc = @(x) x(1,:).^2 + x(2,:).^2;
    %   avgValue = space.average(quadraticFunc);
    %
    % See also:
    %   fem.space.FiniteElementSpace, fem.space.AffineSpace

    properties
        mesh % Geometric discretization object defining the domain
    end

    properties (Dependent)
        nDofs % Total number of degrees of freedom across the entire mesh
    end

    methods
        function obj = MeshSpace(fe, mesh)
            % MESHSPACE Constructor for MeshSpace.
            %
            %   obj = MeshSpace(fe, mesh) creates a finite element function
            %   space defined on the specified geometric mesh using the
            %   given finite element as the local basis.
            %
            % Inputs:
            %   fe - Finite element object defining local basis functions
            %   mesh - Geometric discretization object defining the domain
            %
            % Outputs:
            %   obj - Constructed MeshSpace object

            obj@fem.space.FiniteElementSpace(fe);
            obj.mesh = mesh;
        end

        function n = get.nDofs(obj)
            % GET.NDOFS Get the total number of degrees of freedom.

            n = obj.mesh.nTotalElements * obj.fe.nDofs;
        end

        function Y = evaluate(obj, varargin)
            % EVALUATE Evaluate functions in the mesh space.
            %
            %   Y = evaluate(obj, Z) evaluates using combined input matrix
            %   containing evaluation points and coefficients.
            %
            %   Y = evaluate(obj, X, C) evaluates at physical points X
            %   using coefficient matrix C organized by elements.
            %
            %   Y = evaluate(obj, f, xRef, ...) evaluates function handle f
            %   at reference points on all mesh elements.
            %
            % Inputs:
            %   obj - The MeshSpace object
            %   varargin - Input arguments
            %<   Z - (nDims+nDofs, nPoints) combined matrix [X; C]
            %<   X - (nDims, nPoints) physical evaluation points  
            %<   C - (nDofs, nElements*nFunctions) coefficient matrix
            %<   f - Function handle to evaluate
            %<   xRef - (nDims, nRefPoints) reference points per element
            %
            % Outputs:
            %   Y - Function values at evaluation points

            nArgs = length(varargin);

            if nArgs == 1
                Z = varargin{1};
                fe = obj.fe;
                Y = fe.evaluate(Z);
            elseif nArgs == 2 && ~isa(varargin{1}, 'function_handle')
                [X, C] = varargin{:};
                fe = obj.fe;
                N = obj.mesh.nTotalElements;
                p = size(X, 2);
                m = fe.nDofs;
                C = reshape(C, m, []);
                Y = fe.evaluate(X, C);
                Y = reshape(Y.', N*p, []);
            elseif nArgs >= 2 && isa(varargin{1}, 'function_handle')
                [f, xRef] = varargin{1:2};
                N = obj.mesh.nTotalElements;
                p = size(xRef, 2);
                I = obj.mesh.allElementMultiIndices;
                X = obj.mesh.collocate(xRef, I);
                Y = f(X, varargin{3:end});
                Y = reshape(Y, N*p, []);
            end
        end

        function A = average(obj, f)
            % AVERAGE Compute weighted average of a function over the mesh.
            %
            %   A = average(obj, f) computes the weighted average of the
            %   function f using the integration weights of the finite
            %   element across all mesh elements. This provides a measure
            %   of the function's mean value over the computational domain.
            %
            % Inputs:
            %   obj - The MeshSpace object
            %   f - Function handle to compute average of
            %
            % Outputs:
            %   A - (nElements, nComponents) weighted average values
            
            Th = obj.mesh;
            n = Th.nTotalElements;
            fe = obj.fe;
            D = fe.volumeData;
            w = D.weights;
            I = obj.mesh.allElementMultiIndices;
            xRef = fe.volumeData.nodes;
            X = obj.mesh.collocate(xRef, I);
            A = f(X);
            A = w * reshape(A, length(w), []);
            A = reshape(A, n, []);
        end

        function U = project(obj, f, varargin)
            % PROJECT Project a function onto the finite element space.
            %
            %   U = project(obj, f) projects the function f onto the finite
            %   element space using the appropriate Galerkin projection
            %   method. The projection method (modal or nodal) is
            %   determined by the type of projector associated with the
            %   finite element.
            %
            % Inputs:
            %   obj - The MeshSpace object
            %   f - Function handle to project onto the space
            %   varargin - Input arguments for function handle
            %
            % Outputs:
            %   U - (nDofs, nComponents) projected coefficients

            fe = obj.fe;
            P = fe.projector;
            B = P.basis;

            cls = 'core.function.SeparableFunction';
            core.except.assert( ...
                ~isa(B, cls) || ~iscell(B.factors), ...
                'InvalidBasis', ...
                'All basis functions must have the same type.');

            if isempty(varargin)
                I = obj.mesh.allElementMultiIndices;
                xRef = fe.volumeData.nodes;
                X = obj.mesh.collocate(xRef, I);
                U = f(X);
            else
                U = f(varargin{:});
            end

            if isa(P, 'approx.project.ModalProjector')
                D = fe.volumeData;
                w = D.weights;
                V = D.values;
                F = P.embed(U, V, w);
            end
            
            if isa(P, 'approx.project.NodalProjector')
                F = P.embed(U);
            end

            U = P.project(F);
            U = reshape(U, obj.nDofs, []);
        end
    end
end