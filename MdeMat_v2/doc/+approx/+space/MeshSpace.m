classdef MeshSpace < approx.space.FunctionSpace
    % MESHSPACE Function space defined on a geometric mesh.
    %
    %   MeshSpace handles the mapping between reference elements and
    %   physical elements in the mesh, allowing for efficient computation
    %   of global operations such as L2 projections and weighted averages.
    %   It serves as the foundation for finite element computations on
    %   complex geometries.
    %
    % Examples:
    %   % Create mesh space from element and mesh
    %   space = MeshSpace(element, mesh);
    %   
    %   % Project analytical function onto the mesh space
    %   analyticFunc = @(x) sin(pi*x(1,:)) .* cos(pi*x(2,:));
    %   coefficients = space.project(analyticFunc);
    %   
    %   % Evaluate projected function at specific points
    %   evalPoints = [0.25, 0.75; 0.5, 0.5];
    %   values = space.evaluate([], evalPoints, coefficients);
    %   
    %   % Compute weighted average of a function over the domain
    %   quadraticFunc = @(x) x(1,:).^2 + x(2,:).^2;
    %   avgValue = space.average(quadraticFunc);
    %
    % See also:
    %   approx.space.FunctionSpace, approx.space.AffineSpace

    properties
        mesh % Geometric discretization
    end

    properties (Dependent)
        nTotalElements % Total number of elements
        nLocalDofs % Number of degrees of freedom per element
        nGlobalDofs % Total number of degrees of freedom
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

        function n = get.nTotalElements(obj)
            % GET.NTOTALELEMENTS Get the total number of elements.

            n = obj.mesh.nTotalElements;
        end

        function n = get.nLocalDofs(obj)
            % GET.NLOCALDOFS Get the number of degrees of freedom per
            % element.

            n = obj.element.nDofs;
        end

        function n = get.nGlobalDofs(obj)
            % GET.NGLOBALDOFS Get the total number of degrees of freedom.

            n = obj.nLocalDofs * obj.nTotalElements;
        end
        
        function Y = evaluate(obj, I, varargin)
            % EVALUATE Evaluate functions in the mesh space.
            %
            %   Y = evaluate(obj, Z) evaluates using combined input matrix
            %   containing evaluation points and coefficients.
            %
            %   Y = evaluate(obj, xRef, C) evaluates at reference points
            %   xRef using coefficient matrix C organized by elements.
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

            if isempty(I) && nArgs > 1
                I = obj.mesh.allElementMultiIndices;
            end

            if nArgs == 1
                Z = varargin{1};
                Y = obj.element.evaluate(Z);
            elseif nArgs == 2 && ~isa(varargin{1}, 'function_handle')
                [xRef, C] = varargin{:};
                N = size(I, 1);
                p = size(xRef, 2);
                m = obj.element.nDofs;
                C = reshape(C, m, []);
                Y = obj.element.evaluate(xRef, C);
                Y = reshape(Y.', N*p, []);
            elseif nArgs >= 2 && isa(varargin{1}, 'function_handle')
                [f, xRef] = varargin{1:2};
                N = size(I, 1);
                p = size(xRef, 2);
                X = obj.mesh.collocate(xRef, I);
                Y = f(X, varargin{3:end});
                Y = reshape(Y, N*p, []);
            end
        end

        function A = average(obj, f)
            % AVERAGE Compute weighted average of a function over the mesh.
            %
            %   A = average(obj, f) computes the weighted average of the
            %   function f using the integration weights across all mesh
            %   elements. This provides a measure of the function's mean
            %   value over the computational domain.
            %
            % Inputs:
            %   obj - The MeshSpace object
            %   f - Function handle to compute average
            %
            % Outputs:
            %   A - (nElements, nComponents) weighted average values
            
            n = obj.mesh.nTotalElements;
            w = obj.element.volumeData.weights;
            I = obj.mesh.allElementMultiIndices;
            xRef = obj.element.volumeData.nodes;
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

            P = obj.element.projector;
            B = P.basis;

            cls = 'core.function.SeparableFunction';
            core.except.assert( ...
                ~isa(B, cls) || ~iscell(B.factors), ...
                'InvalidBasis', ...
                'All basis functions must have the same type.');

            if isempty(varargin)
                I = obj.mesh.allElementMultiIndices;
                xRef = obj.element.volumeData.nodes;
                X = obj.mesh.collocate(xRef, I);
                U = f(X);
            else
                U = f(varargin{:});
            end

            if isa(P, 'approx.project.ModalProjector')
                D = obj.element.volumeData;
                w = D.weights;
                V = D.values;
                F = P.embed(U, V, w);
            end
            
            if isa(P, 'approx.project.NodalProjector')
                F = P.embed(U);
            end

            U = P.project(F);
            U = reshape(U, obj.nGlobalDofs, []);
        end
    
    end
end