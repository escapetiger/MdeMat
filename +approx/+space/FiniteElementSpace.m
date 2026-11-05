classdef FiniteElementSpace < approx.space.MeshSpace
    % FINITEELEMENTSPACE Base class for all finite element spaces.
    %
    %   FiniteElementSpace implements the finite-dimensional function space
    %   used in the finite element method, It serves as the foundation for
    %   finite element computations on complex geometries.
    %
    % See also:
    %   approx.space.MeshSpace

    properties
        element % Element object
    end

    properties (Dependent)
        nLocalDofs % Number of local degrees of freedom
        nGlobalDofs % Number of global degrees of freedom
    end

    methods
        function obj = FiniteElementSpace(mesh, element)
            % FINITEELEMENTSPACE Constructor for FiniteElementSpace.
            %
            %   obj = FiniteElementSpace(mesh, element) creates a finite
            %   element space defined by the specified mesh and finite
            %   element.
            %
            % Inputs:
            %   mesh - Mesh object
            %   element - Element object
            %
            % Outputs:
            %   obj - Constructed MeshSpace object

            obj@approx.space.MeshSpace(mesh);
            obj.element = element;
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

        function Y = evaluate(obj, IE, xRef, C)
            % EVALUATE Evaluate finite element functions in the mesh space.
            %
            %   Y = evaluate(obj, IE, xRef, C) evaluates the piecewise
            %   function defined by coefficients @A C at reference points
            %   @a xRef on the specified elements @a IE.
            %
            %   Y = evaluate(obj, [], xRef, C) use all the elements.
            %
            % Inputs:
            %   obj - The MeshSpace object
            %   IE - Element linear indices (ne × 1 or [])
            %   xRef - Reference coordinates (nd × np)
            %   C - Coefficient matrix (nl × ne*nc)
            %
            % Outputs:
            %   Y - Function values at evaluation points (np*ne × nc)

            if isempty(IE)
                IE = (1:obj.nMeshElements).';
            end

            if isempty(xRef)
                xRef = obj.element.volume.nodes;
            end

            ne = size(IE, 1);
            np = size(xRef, 2);
            nl = obj.nLocalDofs;
            C = reshape(C, nl, []);
            Y = obj.element.evaluate(xRef, C); %< ne*nc x np
            Y = reshape(Y.', np*ne, []);
        end

        function Y = differentiate(obj, IE, xRef, C, order)
            % DIFFERENTIATE Differentiate finite element functions in the
            % mesh space.
            %
            %   Y = differentiate(obj, IE, xRef, C, order) differetiates
            %   piecewise function defined by coefficients @A C at
            %   reference points @a xRef on the specified elements @a IE.
            %
            %   Y = differentiate(obj, [], xRef, C, order) use all the
            %   elements.
            %
            % Inputs:
            %   obj - The MeshSpace object
            %   IE - Element linear indices (ne × 1 or [])
            %   xRef - Reference coordinates (nd × np)
            %   C - Coefficient matrix (nl × ne*nc)
            %   order - Derivative order (nd x 1)
            %
            % Outputs:
            %   Y - Function values at evaluation points (np*ne × nc)

            if isempty(IE)
                IE = (1:obj.nMeshElements).';
            end

            if isempty(xRef)
                xRef = obj.element.volume.nodes;
            end

            ne = size(IE, 1);
            np = size(xRef, 2);
            nl = obj.nLocalDofs;
            C = reshape(C, nl, []);
            Y = obj.element.differentiate(xRef, C, order); %< ne*nc x np
            Y = reshape(Y.', np*ne, []);
        end

        function Y = functionEvaluate(obj, IE, xRef, f, varargin)
            % FUNCTIONEVALUATE Evaluate analytical function on mesh
            % elements.
            %
            %   Y = functionEvaluate(obj, IE, xRef, f) evaluates the
            %   analytical function @a f at reference points @a xRef mapped
            %   to physical coordinates on the specified elements @a IE.
            %
            %   Y = functionEvaluate(obj, [], xRef, f, ...) use all
            %   elements.
            %
            %   Y = functionEvaluate(obj, IE, xRef, f, ...) passes
            %   additional arguments to the function @a f.
            %
            % Inputs:
            %   obj - The MeshSpace object
            %   IE - Element linear indices (ne x nd matrix or [])
            %   xRef - Reference coordinates (nd × np matrix)
            %   f - Function handle to evaluate
            %   varargin - Additional arguments for function handle
            %
            % Outputs:
            %   Y - Function values at mapped points (nElements*nRefPoints × nComponents matrix)

            if isempty(IE)
                IE = (1:obj.nMeshElements).';
            end

            core.except.assert(isa(f, 'function_handle'), ...
                'InvalidFunction', 'f must be a function handle.');

            ne = size(IE, 1);
            np = size(xRef, 2);
            X = obj.mesh.collocate(xRef, IE);
            Y = f(X, varargin{:});
            Y = reshape(Y, ne*np, []);
        end

        function A = mean(obj, IE, f, varargin)
            % MEAN Compute weighted mean of a function over specified
            % elements.
            %
            %   A = mean(obj, IE, f) computes the weighted mean of function
            %   @a f using integration weights over the elements specified
            %   by @a IE.
            %
            %   A = mean(obj, [], f) uses all elements.
            %
            %   A = mean(obj, IE, f, ...) passes additional arguments to
            %   function @a f.
            %
            % Inputs:
            %   obj - The MeshSpace object
            %   IE - Element linear indices (ne × nd matrix) or []
            %   f - Function handle to compute mean
            %   varargin - Additional arguments for function handle
            %
            % Outputs:
            %   A - Weighted mean values (ne × nc matrix)

            if isempty(IE)
                IE = (1:obj.nMeshElements).';
            end

            core.except.assert(isa(f, 'function_handle'), ...
                'InvalidFunction', 'f must be a function handle.');

            ne = size(IE, 1);
            np = obj.element.volume.nPoints;
            wRef = obj.element.volume.weights;
            xRef = obj.element.volume.nodes;
            X = obj.mesh.collocate(xRef, IE);
            F = f(X, varargin{:});
            F = reshape(F, np, []);
            A = wRef * F;
            A = reshape(A, ne, []);
        end

        function U = fit(obj, IE, f, varargin)
            % FIT Fit a function within the finite element space.
            %
            %   U = fit(obj, IE, f) fits the function @a f onto the
            %   finite element space over the elements specified by @a IE
            %   using the appropriate fitting method.
            %
            %   U = fit(obj, [], f) uses all elements.
            %
            %   U = fit(obj, IE, f, ...) passes additional arguments to
            %   function @a f.
            %
            % Inputs:
            %   obj - The FiniteElementSpace object
            %   IE - Element linear indices (ne x nd matrix) or []
            %   f - Function handle to fit
            %   varargin - Additional arguments for function handle
            %
            % Outputs:
            %   U - Fitted coefficients (nl*ne × nc matrix)

            if isempty(IE)
                IE = (1:obj.nMeshElements).';
            end

            core.except.assert(isa(f, 'function_handle'), ...
                'InvalidFunction', 'Function f must be a function handle.');

            A = obj.element.approximator;
            B = A.basis;

            cls = 'core.function.SeparableFunction';
            core.except.assert( ...
                ~isa(B, cls) || ~iscell(B.factors), ...
                'InvalidBasis', ...
                'All basis functions must have the same type.');

            ne = size(IE, 1);
            nl = obj.nLocalDofs;
            if isempty(varargin)
                xRef = obj.element.volume.nodes;
                X = obj.mesh.collocate(xRef, IE);
                U = f(X);
            else
                U = f(varargin{:});
            end

            switch A.TYPE
                case 'modal'
                    D = obj.element.volume;
                    F = A.embed(U, D.values, D.weights);
                case 'nodal'
                    F = A.embed(U);
                otherwise
                    core.except.error('UnsupportedProjector', ...
                        'Projector type %s not supported.', class(P));
            end

            U = A.fit(F);
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
        
        function newObj = refine(obj, nLevels)
            % REFINE Create refined finite element space.
            %
            %   newObj = refine(obj, nLevels) creates a refined finite
            %   element space.
            %
            % Inputs:
            %   obj - The FiniteElementSpace object
            %   nLevels - Number of refinement levels (non-negative integer)
            %
            % Outputs:
            %   newObj - Refined FiniteElementSpace object

            newMesh = obj.mesh.refine(nLevels);
            newObj = approx.space.FiniteElementSpace(newMesh, obj.element);
        end
    end
end