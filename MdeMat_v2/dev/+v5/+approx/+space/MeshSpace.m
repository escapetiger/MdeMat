classdef MeshSpace < approx.space.FunctionSpace
    % MESHSPACE Function space defined on a geometric mesh.
    %
    %   MeshSpace handles the mapping between reference elements and
    %   physical elements in the mesh, allowing for efficient computation
    %   of global operations such as L2 projections and weighted averages.
    %   It serves as the foundation for finite element computations on
    %   complex geometries.
    %
    % See also:
    %   approx.space.FunctionSpace, approx.space.AffineSpace

    properties
        mesh % Mesh object
    end

    properties (Dependent)
        nMeshElements % Total number of mesh elements
        nLocalDofs % Number of local degrees of freedom
        nGlobalDofs % Number of global degrees of freedom
    end

    methods
        function obj = MeshSpace(element, mesh)
            % MESHSPACE Constructor for MeshSpace.
            %
            %   obj = MeshSpace(element, mesh) creates a function space
            %   defined on the specified geometric mesh.
            %
            % Inputs:
            %   element - Element object defining local basis functions
            %   mesh - Geometric discretization object defining the domain
            %
            % Outputs:
            %   obj - Constructed MeshSpace object

            obj@approx.space.FunctionSpace(element);
            obj.mesh = mesh;
        end

        function n = get.nMeshElements(obj)
            % GET.NTOTALELEMENTS Get the total number of elements.

            n = obj.mesh.nElements;
        end

        function n = get.nLocalDofs(obj)
            % GET.NLOCALDOFS Get the number of degrees of freedom per
            % element.

            n = obj.element.nDofs;
        end

        function n = get.nGlobalDofs(obj)
            % GET.NGLOBALDOFS Get the total number of degrees of freedom.

            n = obj.nLocalDofs * obj.nMeshElements;
        end

        function Y = evaluate(obj, L, xRef, C)
            % EVALUATE Evaluate finite element functions in the mesh space.
            %
            %   Y = evaluate(obj, L, xRef, C) evaluates the piecewise
            %   function defined by coefficients @A C at reference points
            %   @a xRef on the specified elements @a L.
            %
            %   Y = evaluate(obj, [], xRef, C) use all the elements.
            %
            % Inputs:
            %   obj - The MeshSpace object
            %   L - Element linear indices (ne × 1 vector or [])
            %   xRef - Reference coordinates (nd × np matrix)
            %   C - Coefficient matrix (nl × ne*nc matrix)
            %
            % Outputs:
            %   Y - Function values at evaluation points (np*ne × nc matrix)

            if isempty(L)
                L = (1:obj.nMeshElements).';
            end

            ne = size(L, 1);
            np = size(xRef, 2);
            nl = obj.nLocalDofs;
            C = reshape(C, nl, []);
            Y = obj.element.evaluate(xRef, C); %< ne*nc x np matrix
            Y = reshape(Y.', np*ne, []);
        end

        function Y = functionEvaluate(obj, L, xRef, f, varargin)
            % FUNCTIONEVALUATE Evaluate analytical function on mesh
            % elements.
            %
            %   Y = functionEvaluate(obj, L, xRef, f) evaluates the
            %   analytical function @a f at reference points @a xRef mapped
            %   to physical coordinates on the specified elements @a L.
            %
            %   Y = functionEvaluate(obj, [], xRef, f, ...) use all
            %   elements.
            %
            %   Y = functionEvaluate(obj, L, xRef, f, ...) passes
            %   additional arguments to the function @a f.
            %
            % Inputs:
            %   obj - The MeshSpace object
            %   L - Element linear indices (ne x nd matrix or [])
            %   xRef - Reference coordinates (nd × np matrix)
            %   f - Function handle to evaluate
            %   varargin - Additional arguments for function handle
            %
            % Outputs:
            %   Y - Function values at mapped points (nElements*nRefPoints × nComponents matrix)

            if isempty(L)
                L = (1:obj.nMeshElements).';
            end

            core.except.assert(isa(f, 'function_handle'), ...
                'InvalidFunction', 'f must be a function handle.');

            ne = size(L, 1);
            np = size(xRef, 2);
            X = obj.mesh.collocate(xRef, L);
            Y = f(X, varargin{:});
            Y = reshape(Y, ne*np, []);
        end

        function A = mean(obj, L, f, varargin)
            % MEAN Compute weighted mean of a function over specified
            % elements.
            %
            %   A = mean(obj, L, f) computes the weighted mean of function
            %   @a f using integration weights over the elements specified
            %   by @a L.
            %
            %   A = mean(obj, [], f) uses all elements.
            %
            %   A = mean(obj, L, f, ...) passes additional arguments to
            %   function @a f.
            %
            % Inputs:
            %   obj - The MeshSpace object
            %   L - Element linear indices (ne × nd matrix) or []
            %   f - Function handle to compute mean
            %   varargin - Additional arguments for function handle
            %
            % Outputs:
            %   A - Weighted mean values (ne × nc matrix)

            if isempty(L)
                L = (1:obj.nMeshElements).';
            end

            core.except.assert(isa(f, 'function_handle'), ...
                'InvalidFunction', 'f must be a function handle.');

            ne = size(L, 1);
            np = obj.element.volumeData.nPoints;
            wRef = obj.element.volumeData.weights;
            xRef = obj.element.volumeData.nodes;
            X = obj.mesh.collocate(xRef, L);
            F = f(X, varargin{:});
            F = reshape(F, np, []);
            A = wRef * F;
            A = reshape(A, ne, []);
        end

        function U = project(obj, L, f, varargin)
            % PROJECT Project a function onto the finite element space.
            %
            %   U = project(obj, L, f) projects the function @a f onto the
            %   finite element space over the elements specified by @a L
            %   using the appropriate Galerkin projection method.
            %
            %   U = project(obj, [], f) uses all elements.
            %
            %   U = project(obj, L, f, ...) passes additional arguments to
            %   function @a f.
            %
            % Inputs:
            %   obj - The MeshSpace object
            %   L - Element linear indices (ne x nd matrix) or []
            %   f - Function handle to project onto the space
            %   varargin - Additional arguments for function handle
            %
            % Outputs:
            %   U - Projected coefficients (nl*ne × nc matrix)

            if isempty(L)
                L = (1:obj.nMeshElements).';
            end

            core.except.assert(isa(f, 'function_handle'), ...
                'InvalidFunction', 'Function f must be a function handle.');

            P = obj.element.projector;
            B = P.basis;

            cls = 'core.function.SeparableFunction';
            core.except.assert( ...
                ~isa(B, cls) || ~iscell(B.factors), ...
                'InvalidBasis', ...
                'All basis functions must have the same type.');

            ne = size(L, 1);
            nl = obj.nLocalDofs;
            if isempty(varargin)
                xRef = obj.element.volumeData.nodes;
                X = obj.mesh.collocate(xRef, L);
                U = f(X);
            else
                U = f(varargin{:});
            end

            switch P.TYPE
                case 'modal'
                    D = obj.element.volumeData;
                    F = P.embed(U, D.values, D.weights);
                case 'nodal'
                    F = P.embed(U);
                otherwise
                    core.except.error('UnsupportedProjector', ...
                        'Projector type %s not supported.', class(P));
            end

            U = P.project(F);
            U = reshape(U, nl*ne, []);
        end

        function G = localToGlobal(obj, E, L)
            % LOCALTOGLOBAL Map local DoF indices to global DoF indices.
            %
            %   G = localToGlobal(obj, E, L) maps local DoF indices L to
            %   global DoF indices for the specified elements E.
            %
            % Inputs:
            %   obj - The DofMap object
            %   E - Element linear indices (ne x 1 vector)
            %   L - Local DoF indices (nl x 1 vector or [])
            %
            % Outputs:
            %   G - Global DoF indices

            nl = obj.nLocalDofs;

            if isempty(E)
                E = (1:obj.nMeshElements);
            end

            if isempty(L)
                L = (1:nl);
            end
            G = L(:) + (E(:).' - 1) * nl;
            G = G(:);
        end

        function [E, L] = globalToLocal(obj, G)
            % GLOBALTOLOCAL Map global DoF indices to local DoF indices.
            %
            %   [E, L] = globalToLocal(obj, G) maps global DoF indices @a G
            %   to local DoF indices.
            %
            % Inputs:
            %   obj - The DofMap object
            %   G - Global DoF indices (ng x 1 vector)
            %
            % Outputs:
            %   E - Element linear indices (ng x 1 vector)
            %   L - Local DoF indices (ng x 1 vector)

            nl = obj.nLocalDofs;
            E = ceil(G/nl);
            L = G - (E - 1) * nl;
        end
    end
end