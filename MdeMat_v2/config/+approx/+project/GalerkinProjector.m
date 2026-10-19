classdef GalerkinProjector < handle
    % GALERKINPROJECTOR Base class for all Galerkin-type projectors.
    %
    %   GalerkinProjector provides the foundation for Galerkin
    %   approximation methods where functions are projected onto
    %   finite-dimensional spaces using mass operators to solve linear
    %   systems. The class implements the core Galerkin methodology for
    %   function approximation.
    %
    %   The Galerkin method seeks an approximation \f$u_h \in V_h\f$ of a
    %   function \f$u\f$ such that
    %   \f[
    %   \langle u_h, v_h\rangle = \langle f, v_h\rangle, \forall v_h ∈ V_h,
    %   \f]
    %   where \f$V_h\f$ is the finite-dimensional function space spanned
    %   by the basis.
    %
    %   This class defines the interface for concrete implementations such
    %   as modal and nodal projectors, which differ in how they construct
    %   the mass operator and embed data into the latent space.
    %
    % See Also:
    %   approx.project.ModalProjector, approx.project.NodalProjector

    properties
        basis % Basis functions that span the approximation space
        mass  % Mass operator matrix for the Galerkin system
    end

    properties (Dependent)
        nDofs % Number of degrees of freedom
        isWellPosed % Whether the Galerkin system is well-posed
        isOverdetermined % Whether the Galerkin system is overdetermined
        isUnderdetermined % Whether the Galerkin system is underdetermined
    end

    methods
        function obj = GalerkinProjector(basis)
            % GALERKINPROJECTOR Constructor for GalerkinProjector.
            %
            %   obj = GalerkinProjector(basis) creates a new
            %   GalerkinProjector object with the specified basis
            %   functions. The mass operator is initialized as empty and
            %   must be set by calling setMass().
            %
            % Inputs:
            %   basis - Basis functions that span the approximation space
            %
            % Outputs:
            %   obj - The constructed GalerkinProjector object
            
            core.except.assert(~isempty(basis), 'InvalidInput', ...
                'Basis functions cannot be empty.');
            
            obj.basis = basis;
            obj.mass = [];
        end        

        function n = get.nDofs(obj)
            % GET.NDOFS Returns the number of degrees of freedom.

            n = obj.basis.nCodims;
        end

        function Y = evaluate(obj, varargin)
            % EVALUATE Evaluate functions in the approximation space at
            % given points.
            %
            %   Y = evaluate(obj, Z) evaluates functions represented by
            %   coefficients embedded in the augmented vector Z = [X; C]
            %   at evaluation points X.
            %
            %   Y = evaluate(obj, X, C) evaluates functions with
            %   coefficient matrix C at evaluation points X.
            %
            % Inputs:
            %   obj - The GalerkinProjector object
            %   varargin - Input arguments
            %<   Z - (n+m, p) matrix representing [X; C]
            %<   X - (n, p) matrix representing evaluation points
            %<   C - (m, q) matrix representing coefficients
            %
            % Outputs:
            %   Y - Function values ((1 x p) or (q x p))
            
            switch length(varargin)
                case 1
                    n = obj.basis.nDims;
                    Z = varargin{1};
                    Z = reshape(Z, size(Z, 1), []);
                    
                    core.except.assert(size(Z, 1) > n, 'InvalidInput', ...
                        'Input Z must have more rows than basis dimensions.');
                    
                    C = Z((n+1):end, :);
                    X = Z(1:n, :);
                    B = obj.basis.evaluate(X);
                    Y = dot(B, C, 1);
                    
                case 2
                    [X, C] = varargin{:};
                    X = reshape(X, size(X, 1), []);
                    C = reshape(C, size(C, 1), []);
                    B = obj.basis.evaluate(X);
                    Y = C.' * B;
                    
                otherwise
                    core.except.assert(0, 'InvalidInput', ...
                        'Invalid number of input arguments.');
            end
        end

        function U = project(obj, F)
            % PROJECT Project data onto the function space.
            %
            %   U = project(obj, F) projects data from the latent space
            %   onto the function space by solving the Galerkin system M*U
            %   = F, where M is the mass operator. Handles well-posed and
            %   overdetermined systems automatically.
            %
            % Inputs:
            %   obj - The GalerkinProjector object
            %   F - Data in the latent space (m x ...)
            %
            % Outputs:
            %   U - Coefficients of the basis functions (m x ...)
            
            M = obj.mass;
            core.except.assert(~isempty(M), 'EmptyMass', ...
                'The mass operator has not been set. Call setMass() first.');

            m = size(M, 1);
            core.except.assert(size(F, 1) == m, 'InvalidInput', ...
                'Data dimension (%d) must match mass operator rows (%d).', ...
                size(F, 1), m);
                        
            s = size(F);
            F = reshape(F, m, []);
            
            if obj.isWellPosed
                U = M \ F;
            elseif obj.isOverdetermined
                U = (M.' * M) \ (M.' * F);
            else
                core.except.assert(0, 'UnderdeterminedSystem', ...
                    'Mass operator cannot be underdetermined.');
            end
            
            U = reshape(U, [size(U, 1), s(2:end)]);
        end
        
        function TF = get.isWellPosed(obj)
            % GET.ISWELLPOSED Check if the Galerkin system is well-posed.
            
            M = obj.mass;
            core.except.assert(~isempty(M), 'EmptyMass', ...
                'The mass operator has not been set.');
            
            TF = size(M, 1) == size(M, 2);
        end
        
        function TF = get.isOverdetermined(obj)
            % GET.ISOVERDETERMINED Check if the Galerkin system is
            % overdetermined.
            
            M = obj.mass;
            core.except.assert(~isempty(M), 'EmptyMass', ...
                'The mass operator has not been set.');
            
            TF = size(M, 1) > size(M, 2);
        end
        
        function TF = get.isUnderdetermined(obj)
            % ISUNDERDETERMINED Check if the Galerkin system is
            % underdetermined.
            
            M = obj.mass;
            core.except.assert(~isempty(M), 'EmptyMass', ...
                'The mass operator has not been set.');
            
            TF = size(M, 1) < size(M, 2);
        end
    end
    
    methods (Abstract)
        % SETMASS Set up the mass operator for the Galerkin method.
        obj = setMass(obj, varargin)

        % EMBED Embed data into the latent space.
        F = embed(obj, varargin)
    end
end