classdef DiffusionGridAssembly < fem.assembly.GridAssembly
    % DIFFUSIONGRIDASSEMBLY Grid-based assembly for diffusion operators.
    %
    %   DiffusionGridAssembly provides assembly methods for diffusion-type
    %   operators on structured grid meshes. This class combines auxiliary
    %   and primal variable assemblys to construct complete diffusion
    %   operators with various jump treatments and boundary conditions.
    %
    %   The assembly supports both implicit and explicit jump treatments
    %   for handling discontinuities at element interfaces, which is
    %   essential for discontinuous Galerkin and mixed finite element
    %   methods.
    %
    % Examples:
    %   % Create diffusion assembly with implicit jumps
    %   assembly = DiffusionGridAssembly(fe, mesh, op, 1, 2, 1, 1);
    %   
    %   % Assemble isotropic diffusion operator
    %   diffusionMatrix = assembly.hessian(eye(2));
    %
    %   % Assemble anisotropic diffusion with boundary conditions
    %   anisotropyTensor = [2.0, 0.1; 0.1, 0.5];
    %   bcVector = assembly.hessianBc(anisotropyTensor, @bcFunc, @jumpFunc);
    %
    % See also:
    %   fem.assembly.GridAssembly, fem.assembly.AuxiliaryGridAssembly,
    %   fem.assembly.PrimalGridAssembly

    properties
        aux      % Auxiliary variable assembly structure
        prm      % Primal variable assembly structure  
        jumpType % Jump treatment type (1=implicit, 2=explicit)
    end

    methods
        function obj = DiffusionGridAssembly(fe, mesh, op, bcType, auxFlux, prmFlux, jumpType)
            % DIFFUSIONGRIDASSEMBLY Constructor for DiffusionGridAssembly.
            %
            %   obj = DiffusionGridAssembly(fe, mesh, op, bcType, auxFlux,
            %   prmFlux, jumpType) creates an assembly for diffusion
            %   operators using mixed formulations with auxiliary and
            %   primal variables.
            %
            % Inputs:
            %   fe - Finite element object
            %   mesh - Grid mesh object (must be approx.mesh.Grid)
            %   op - Finite element operator object
            %   bcType - Boundary condition type (0=periodic, 1=Dirichlet)
            %   auxFlux - Flux type for auxiliary variables (1=left, 2=right, 3=central)
            %   prmFlux - Flux type for primal variables (1=left, 2=right, 3=central)
            %   jumpType - Jump treatment (1=implicit, 2=explicit)
            %
            % Outputs:
            %   obj - Constructed DiffusionGridAssembly object

            core.except.assert(isa(mesh, 'approx.mesh.Grid'), ...
                'InvalidInput', ...
                'Finite element space must be based on grid.');

            obj@fem.assembly.GridAssembly(fe, mesh, op, bcType);
            obj.aux.asm = fem.assembly.AuxiliaryGridAssembly(fe, mesh, op, bcType, auxFlux);
            obj.aux.mat = cell(1, obj.nDims);
            obj.prm.asm = fem.assembly.PrimalGridAssembly(fe, mesh, op, bcType, prmFlux);
            obj.prm.mat = cell(1, obj.nDims);
            obj.jumpType = jumpType;
        end

        function A = hessian(obj, c)
            % HESSIAN Assemble diffusion Hessian operator.
            %
            %   A = hessian(obj, c) assembles the complete diffusion
            %   operator matrix using the specified diffusion tensor. The
            %   method combines auxiliary and primal variable contributions
            %   with appropriate jump terms based on the jump type.
            %
            % Inputs:
            %   obj - The DiffusionGridAssembly object
            %   c - Diffusion tensor (nDims×nDims matrix)
            %
            % Outputs:
            %   A - Assembled diffusion operator sparse matrix
            %
            % Examples:
            %   % Isotropic diffusion
            %   isotropicMatrix = assembly.hessian(eye(2));
            %
            %   % Anisotropic diffusion
            %   tensor = [2.0, 0.5; 0.5, 1.0];
            %   anisotropicMatrix = assembly.hessian(tensor);

            for i = 1:obj.nDims
                obj.aux.mat{i} = obj.aux.asm.partial(i, 1);
            end
            for i = 1:obj.nDims
                obj.prm.mat{i} = obj.prm.asm.partial(i, 1);
            end
            n = obj.nGlobalDofs;
            A = sparse(n, n);
            for i = 1:obj.nDims
                for j = 1:obj.nDims
                    if c(i, j) ~= 0
                        A = A + c(i, j) * obj.aux.mat{i} * obj.prm.mat{j};
                        if obj.bcType && obj.jumpType == 1
                            B = obj.aux.asm.outerJump(j, c(i, j));
                            A = A + c(i, j) * B;
                        end
                    end
                end
            end
        end

        function b = hessianBc(obj, c, f, g, varargin)
            % HESSIANBC Assemble boundary condition vector for diffusion
            % operator.
            %
            %   b = hessianBc(obj, c, f, g) assembles the boundary
            %   condition contribution vector for the diffusion operator
            %   with specified boundary functions and diffusion tensor.
            %
            % Inputs:
            %   obj - The DiffusionGridAssembly object
            %   c - Diffusion tensor (nDims×nDims matrix)
            %   f - Boundary condition function handle
            %   g - Extraction function handle (optional for explicit jumps)
            %   varargin - Additional arguments for boundary functions
            %
            % Outputs:
            %   b - Assembled boundary condition sparse vector
            %
            % Examples:
            %   % Simple Dirichlet boundary conditions
            %   bcVector = assembly.hessianBc(eye(2), @(x) sin(x), []);
            %
            %   % With explicit jump treatment
            %   bcVector = assembly.hessianBc(tensor, @bcFunc, @jumpFunc, param1);

            n = obj.nGlobalDofs;
            b = sparse(n, 1);
            for i = 1:obj.nDims
                for j = 1:obj.nDims
                    b = b + obj.partialBc(i, j, c(i, j), f, g, varargin{:});
                end
            end
        end

        function b = partialBc(obj, i, j, c, f, g, varargin)
            % PARTIALBC Assemble partial boundary condition contribution.
            %
            %   b = partialBc(obj, i, j, c, f, g) assembles the boundary
            %   condition contribution for a specific component of the
            %   diffusion tensor, handling both Dirichlet conditions and
            %   jump terms as appropriate.
            %
            % Inputs:
            %   obj - The DiffusionGridAssembly object
            %   i - First tensor index (row)
            %   j - Second tensor index (column)
            %   c - Tensor component value
            %   f - Boundary condition function handle
            %   g - Jump function handle (for explicit jumps)
            %   varargin - Additional arguments for boundary functions
            %
            % Outputs:
            %   b - Boundary condition contribution sparse vector

            ne = obj.nTotalElements;

            %< ------------------------------------------------------------
            %< COEFFICIENTS
            %< ------------------------------------------------------------
            if isa(obj.mesh, 'approx.mesh.UniformGrid')
                h = obj.mesh.spacings{j};
                c = ones(1, ne) .* c ./ h;
            end

            if isa(obj.mesh, 'approx.mesh.NonuniformGrid')
                h = obj.mesh.spacings;
                [h{:}] = ndgrid(h{:});
                c = c ./ h{j}(:);
            end

            c = c(:);

            %< ------------------------------------------------------------
            %< FLUX CONTRIBUTION
            %< ------------------------------------------------------------
            if isa(obj.fe, 'fem.element.L2FiniteElement')
                T1 = obj.assembleLeftDirichletBoundaryConditionTerms(j, c, f, varargin{:});
                T2 = obj.assembleRightDirichletBoundaryConditionTerms(j, c, f, varargin{:});
            end

            %< ------------------------------------------------------------
            %< TRIPLETS TO SPARSE MATRIX
            %< ------------------------------------------------------------
            n = obj.nGlobalDofs;
            T = [T1; T2];

            b = obj.sparseFromTriplets(T, n, 1);
            b = obj.aux.mat{i} * b;
            if obj.jumpType == 1
                if obj.aux.asm.fluxType == 1
                    b = b + obj.sparseFromTriplets(T1, n, 1);
                elseif obj.aux.asm.fluxType == 2
                    b = b + obj.sparseFromTriplets(T2, n, 1);
                elseif obj.aux.asm.fluxType == 3
                    b = b + b / 2;
                end
            end
            if ~isempty(g) && obj.jumpType == 2
                if obj.prm.asm.fluxType == 1
                    TT = obj.assembleLeftDirichletBoundaryJumpTerms(j, c, f, g, varargin{:});
                elseif obj.prm.asm.fluxType == 2
                    TT = obj.assembleRightDirichletBoundaryJumpTerms(j, c, f, g, varargin{:});
                elseif obj.prm.asm.fluxType == 3
                    TT1 = obj.assembleLeftDirichletBoundaryJumpTerms(j, c/2, f, g, varargin{:});
                    TT2 = obj.assembleRightDirichletBoundaryJumpTerms(j, c/2, f, g, varargin{:});
                    TT = [TT1; TT2];
                end
                b = b + obj.sparseFromTriplets(TT, n, 1);
            end
        end
    end
end