classdef EllipticAssembly < mdemat.common.fem.AdjointAssembly
    % FINITEELEMENTELLIPTICASSEMBLY Assembly for finite element elliptic
    % operators.
    
    properties
        Cache % Cached data
    end
    
    methods
        function obj = EllipticAssembly(space, bcType, penaltyType)
            % FINITEELEMENTELLIPTICASSEMBLY Constructor for elliptic
            % assembly.
            %
            %   obj = EllipticAssembly(space, bcType,
            %   penaltyType) creates an elliptic assembly object with
            %   specified finite element @a space, boundary condition @a
            %   bcType, and penalty formulation @a penaltyType. Inherits
            %   from AdjointAssembly and initializes the cache
            %   for efficient matrix-vector operations.
            
            arguments
                space mdemat.common.fem.FiniteElementSpace
                bcType{mustBeMember(bcType, {'periodic', 'dirichlet'})}
                penaltyType{mustBePenaltyType(penaltyType)}
            end
            
            obj@mdemat.common.fem.AdjointAssembly(space, bcType, penaltyType);
            obj.reset();
        end
        
        function obj = reset(obj)
            % RESET Clear the assembly cache.
            %
            %   obj = reset(obj) clears the internal cache structure used
            %   to store intermediate assembly results for efficient
            %   computation. This method should be called when the finite
            %   element space or boundary conditions change to ensure cache
            %   consistency.
            
            arguments
                obj mdemat.common.fem.EllipticAssembly
            end
            
            obj.Cache = struct();
        end
        
        function [A, b] = forceCompletion(obj, A, b)
            % FORCECOMPLETION Force completion of singular elliptic system.
            %
            %   [A, b] = forceCompletion(obj, A, b) augments the matrix @a
            %   A and vector @a b to resolve singular elliptic systems that
            %   arise with periodic boundary conditions. Adds a
            %   constraint that enforces zero mean of the solution by
            %   appending a Lagrange multiplier equation to the system.
            %
            %   The augmented system has the form:
            %   \f[
            %       \begin{bmatrix} A & c \\ \ell^T & 0 \end{bmatrix}
            %       \begin{bmatrix} u \\ \lambda \end{bmatrix} =
            %       \begin{bmatrix} b \\ 0 \end{bmatrix}
            %   \f]
            %
            %   Here c is the coefficient representation of the constant
            %   source. The physical integration weights ell enforce zero
            %   mean, ell^T*u=0, and include each element Jacobian. This is
            %   not a Hermitian transpose: the integral is linear in u.
            
            arguments
                obj mdemat.common.fem.EllipticAssembly
                A{mustBeNumeric}
                b{mustBeNumeric}
            end
            
            source = mdemat.common.fem.SourceAssembly(obj.Space);
            c = source.assembleVector(1);
            integral = obj.Space.integrationWeights();
            A = [A, c; integral.', 0];
            b = [b; zeros(1, size(b, 2))];
        end
        
        function G = assembleMatrix(obj, field)
            % ASSEMBLEMATRIX Assemble elliptic operator matrix.
            %
            %   G = assembleMatrix(obj, field) assembles the global matrix
            %   for the diffusion generator \f$ \nabla \cdot (D \nabla u)
            %   \f$ where @a field represents the diffusion tensor \f$ D
            %   \f$. The assembly combines volume and flux contributions
            %   using the discontinuous Galerkin method with interior
            %   penalty stabilization.
            %
            %   For Dirichlet boundary conditions, intermediate matrices
            %   are cached in obj.Cache for efficient vector assembly.
            %
            %   The matrix is the mass-normalized diffusion generator.
            %   Its sign is opposite to the positive elliptic stiffness:
            %
            %   \f[
            %       (G u,v) = -\sum_K \int_K D \nabla u \cdot \nabla v \, dx
            %                + \text{DG interface and boundary terms}.
            %   \f]
            
            arguments
                obj mdemat.common.fem.EllipticAssembly
                field{mustBeNumeric}
            end
            
            nd = obj.Space.NDims;
            
            field = reshape(field, nd, nd, []);
            
            GV = obj.assembleConstantVolumeMatrix(field);
            GF = obj.assembleConstantFluxMatrix(field);
            G = obj.createZeros();
            if strcmpi(obj.BcType, 'dirichlet')
                obj.Cache.G12 = cell(1, nd);
            end
            for d = 1:nd
                G12 = GF{1, 2}{d} - GV{1, 2}{d};
                G21 = GF{2, 1}{d} - GV{2, 1}{d};
                G11 = GF{1, 1}{d};
                % G11 is the negative jump generator, so it adds
                % dissipation to diffusion and changes sign in -Delta.
                G = G + G12 * G21 + G11;
                if strcmpi(obj.BcType, 'dirichlet')
                    obj.Cache.G12{d} = G12;
                end
            end
        end
        
        function G = assembleVector(obj, field, func, options)
            % ASSEMBLEVECTOR Assemble vector contribution for elliptic
            % operator.
            %
            %   G = assembleVector(obj, field, func) assembles the vector
            %   contribution for the elliptic operator with diffusion @a
            %   field and boundary function @a func. Returns the assembled
            %   right-hand side vector.
            %
            %   G = assembleVector(obj, field, func, options) allows
            %   additional arguments to be passed to the function
            %   evaluation through @a options.args.
            
            arguments
                obj mdemat.common.fem.EllipticAssembly
                field{mustBeNumeric}
                func{mustBeA(func, 'function_handle')}
                options.mode = ''
                options.args = {}
            end
            
            if strcmpi(obj.BcType, 'periodic')
                G = obj.createZeros(nCols=1);
                return;
            end
            
            nd = obj.Space.NDims;
            
            field = reshape(field, nd, nd, []);
            
            GF = obj.assembleConstantFluxVector(field, func, options.mode, options.args{:});
            
            G = obj.createZeros(nCols=1);
            for d = 1:nd
                G = G + obj.Cache.G12{d} * GF{2, 1}{d} + GF{1, 1}{d};
            end
        end
    end
end
