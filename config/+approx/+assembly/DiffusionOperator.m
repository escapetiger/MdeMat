classdef DiffusionOperator < approx.assembly.Assembly
    % DIFFUSIONOPERATOR Assembly for diffusion operators in LDG methods.
    %
    %   DiffusionOperator provides high-level assembly of diffusion
    %   operators using discontinuous Galerkin discretizations. It
    %   implements the Local Discontinuous Galerkin (LDG) method by
    %   combining auxiliary and primal flux formulations to handle
    %   second-order elliptic operators.
    %
    %   The operator discretizes terms of the form:
    %
    %   \f[
    %     -\nabla \cdot (D \nabla u)
    %   \f]
    %
    %   where \f$D\f$ is a diffusion tensor and \f$u\f$ is the diffused
    %   quantity.
    %
    %   The LDG method introduces auxiliary variables for gradients and
    %   uses both auxiliary and primal flux treatments. Jump terms handle
    %   the enforcement of boundary conditions, with options for implicit
    %   or explicit treatment of Dirichlet conditions.
    %
    % Examples:
    %   % Create diffusion operator with Dirichlet boundaries
    %   space = approx.space.MeshSpace(mesh, element);
    %   operator = approx.element.L2ElementOperator(element);
    %   diffOp = DiffusionOperator(space, operator, 'dirichlet', ...
    %                             'central', 'central', 'implicit');
    %
    %   % Assemble diffusion matrix for isotropic diffusion
    %   diffusionTensor = eye(2);  % Isotropic diffusion in 2D
    %   A = diffOp.linear(diffusionTensor);
    %
    %   % Assemble boundary condition terms
    %   boundaryFunction = @(x) x(1,:).^2 + x(2,:).^2;
    %   interiorFunction = @(x) 2*x(1,:) + 2*x(2,:);  % Gradient extraction
    %   B = diffOp.linearBc(diag(diffusionTensor), boundaryFunction, interiorFunction);
    %
    % See also:
    %   approx.assembly.Assembly, approx.assembly.VolumeAssembly,
    %   approx.assembly.FluxAssembly, approx.assembly.ConvectionOperator

    properties
        volume % Volume assembly for interior terms (VolumeAssembly)
        auxFlux % Auxiliary flux assembly (FluxAssembly)
        prmFlux % Primal flux assembly (FluxAssembly)
        auxMat % Auxiliary operator matrices (cell array)
        prmMat % Primal operator matrices (cell array)
        jmpMat % Jump operator matrices (cell array)
        trMat % Trace matrices (cell array)
        trJmpMat % Trace jump matrices (cell array)
        jumpType % Jump treatment type: {'implicit', 'explicit'}
    end

    properties (Dependent)
        hasDirichletBoundaryJump % True if using any Dirichlet jump treatment
        isImplicitDirichletBoundaryJump % True if using implicit Dirichlet jump treatment
        isExplicitDirichletBoundaryJump % True if using explicit Dirichlet jump treatment
    end

    methods
        function obj = DiffusionOperator(space, operator, bcType, auxFluxType, prmFluxType, jumpType)
            % DIFFUSIONOPERATOR Constructor for DiffusionOperator.
            %
            %   obj = DiffusionOperator(space, operator, bcType,
            %   auxFluxType, prmFluxType, jumpType) creates a diffusion
            %   operator with the specified components for Local
            %   Discontinuous Galerkin discretization.
            %
            % Inputs:
            %   space - MeshSpace object
            %   operator - ElementOperator object
            %   bcType - Boundary condition type: {'periodic', 'dirichlet'}
            %   auxFluxType - Auxiliary flux type: {'left', 'right', 'central'}
            %   prmFluxType - Primal flux type: {'left', 'right', 'central'}
            %   jumpType - Jump treatment: {'implicit', 'explicit'}
            %
            % Outputs:
            %   obj - Constructed DiffusionOperator object

            obj@approx.assembly.Assembly(space, operator);
            obj.volume = approx.assembly.VolumeAssembly(space, operator);
            obj.auxFlux = approx.assembly.FluxAssembly(space, operator, bcType, auxFluxType);
            obj.prmFlux = approx.assembly.FluxAssembly(space, operator, bcType, prmFluxType);
            if strcmpi(bcType, 'dirichlet')
                obj.auxFlux.setFluxMode('auxiliary');
                obj.prmFlux.setFluxMode('primal');
            end
            obj.auxMat = cell(1, obj.space.nDims);
            obj.prmMat = cell(1, obj.space.nDims);
            obj.jmpMat = cell(1, obj.space.nDims);
            obj.trMat = cell(1, 2*obj.space.nDims);
            obj.trJmpMat = cell(1, 2*obj.space.nDims);
            obj.jumpType = jumpType;
        end

        function TF = get.hasDirichletBoundaryJump(obj)
            % GET.HASDIRICHLETBOUNDARYJUMP Check for Dirichlet jump
            % treatment.
            TF = obj.isImplicitDirichletBoundaryJump || obj.isExplicitDirichletBoundaryJump;
        end

        function TF = get.isImplicitDirichletBoundaryJump(obj)
            % GET.ISIMPLICITDIRICHLETBOUNDARYJUMP Check for implicit
            % Dirichlet jump treatment.

            TF = strcmpi(obj.auxFlux.bcType, 'dirichlet') && strcmpi(obj.jumpType, 'implicit');
        end

        function TF = get.isExplicitDirichletBoundaryJump(obj)
            % GET.ISEXPLICITDIRICHLETBOUNDARYJUMP Check for explicit
            % Dirichlet jump treatment.

            TF = strcmpi(obj.auxFlux.bcType, 'dirichlet') && strcmpi(obj.jumpType, 'explicit');
        end

        function A = linear(obj, coe)
            % LINEAR Assemble linear diffusion operator matrix.
            %
            %   A = linear(obj, coe) assembles the diffusion operator
            %   matrix using the Local Discontinuous Galerkin method. The
            %   coefficient matrix @a coe represents the diffusion tensor.
            %
            % Inputs:
            %   obj - The DiffusionOperator object
            %   coe - Diffusion tensor (nDims x nDims matrix)
            %
            % Outputs:
            %   A - Assembled diffusion operator matrix (sparse)

            n = obj.space.nGlobalDofs;
            h = obj.space.mesh.spacings;

            for dim = 1:obj.space.nDims
                fa = isempty(obj.auxMat{dim});
                fp = isempty(obj.prmMat{dim});
                fj = isempty(obj.jmpMat{dim}) && obj.isImplicitDirichletBoundaryJump;
                if fa || fp || fj
                    C = obj.volume.scaleConstant(dim, 1);
                    T1 = obj.volume.assembleVolumePartial(dim, C);
                end

                if fa
                    T2 = obj.auxFlux.assembleFluxPartial(dim, C);
                    T = [T1; T2];
                    obj.auxMat{dim} = core.linalg.sparseFromTriplets(T, n, n);
                end

                if fp
                    T2 = obj.prmFlux.assembleFluxPartial(dim, C);
                    T = [T1; T2];
                    obj.prmMat{dim} = core.linalg.sparseFromTriplets(T, n, n);
                end

                if fj
                    cL = -1 ./ h{dim};
                    cR = 1 ./ h{dim};
                    switch lower(obj.auxFlux.fluxType)
                        case 'left'
                            T = obj.auxFlux.boundary.assembleImplicitLeftJump(dim, cL.*C);
                        case 'right'
                            T = obj.auxFlux.boundary.assembleImplicitRightJump(dim, -cR.*C);
                        case 'central'
                            T1 = obj.auxFlux.boundary.assembleImplicitLeftJump(dim, cL.*C/2);
                            T2 = obj.auxFlux.boundary.assembleImplicitRightJump(dim, -cR.*C/2);
                            T = [T1; T2];
                    end
                    obj.jmpMat{dim} = core.linalg.sparseFromTriplets(T, n, n);
                end
            end

            A = sparse(n, n);
            for i = 1:obj.space.nDims
                for j = 1:obj.space.nDims
                    if coe(i, j) == 0, continue; end
                    A = A + coe(i, j) * obj.auxMat{i} * obj.prmMat{j};
                end
            end

            if obj.isImplicitDirichletBoundaryJump
                for i = 1:obj.space.nDims
                    for j = 1:obj.space.nDims
                        if coe(i, j) == 0, continue; end
                        A = A + coe(i, j) * obj.jmpMat{j};
                    end
                end
            end
        end

        function A = linearBc(obj, coe, f, g, varargin)
            % LINEARBC Assemble boundary condition terms for linear
            % diffusion.
            %
            %   A = linearBc(obj, coe, f, g, varargin) assembles the
            %   boundary condition contribution matrix for linear diffusion
            %   with Dirichlet boundary conditions. The treatment depends
            %   on whether implicit or explicit jump handling is used.
            %
            % Inputs:
            %   obj - The DiffusionOperator object
            %   coe - Diagonal diffusion coefficients (nDims x 1)
            %   f - Boundary condition function handle
            %   g - Interior extraction function handle (for explicit jumps)
            %   varargin - Input arguments
            %
            % Outputs:
            %   A - Boundary condition contribution matrix (sparse)

            n = obj.space.nGlobalDofs;
            h = obj.space.mesh.spacings;

            for dim = 1:obj.space.nDims
                C = obj.volume.scaleConstant(dim, 1);
                T = obj.auxFlux.boundary.assembleTrace(2*dim-1, C, f, varargin{:});
                obj.trMat{2*dim-1} = core.linalg.sparseFromTriplets(T, n, max(T(:, 2)));
                T = obj.auxFlux.boundary.assembleTrace(2*dim, C, f, varargin{:});
                obj.trMat{2*dim} = core.linalg.sparseFromTriplets(T, n, max(T(:, 2)));
            end

            A = 0;
            for i = 1:obj.space.nDims
                for j = 1:obj.space.nDims
                    if coe(i, j) == 0, continue; end
                    A = A + coe(i, j) * obj.auxMat{i} * (obj.trMat{2*j} - obj.trMat{2*j-1});
                end
            end

            if obj.isImplicitDirichletBoundaryJump
                for dim = 1:obj.space.nDims
                    cL = -1 ./ h{dim};
                    cR = 1 ./ h{dim};
                    C = obj.volume.scaleConstant(dim, 1);
                    T = obj.auxFlux.boundary.assembleTrace(2*dim-1, -cL.*C, f, varargin{:});
                    obj.trJmpMat{2*dim-1} = core.linalg.sparseFromTriplets(T, n, max(T(:, 2)));
                    T = obj.auxFlux.boundary.assembleTrace(2*dim, cR.*C, f, varargin{:});
                    obj.trJmpMat{2*dim} = core.linalg.sparseFromTriplets(T, n, max(T(:, 2)));
                end

                for i = 1:obj.space.nDims
                    for j = 1:obj.space.nDims
                        if coe(i, j) == 0, continue; end
                        switch obj.auxFlux.fluxType
                            case 'left'
                                J = obj.trJmpMat{2*j-1};
                            case 'right'
                                J = obj.trJmpMat{2*j};
                            case 'central'
                                J = (obj.trJmpMat{2*j-1} + obj.trJmpMat{2*j}) / 2;
                        end
                        A = A + coe(i, j) * J;
                    end
                end
            elseif obj.isExplicitDirichletBoundaryJump
                for dim = 1:obj.space.nDims
                    cL = -1 ./ h{dim};
                    cR = 1 ./ h{dim};
                    C = obj.volume.scaleConstant(dim, 1);
                    T = obj.auxFlux.boundary.assembleJump(2*dim-1, cL*C, f, g, varargin{:});
                    obj.trJmpMat{2*dim-1} = core.linalg.sparseFromTriplets(T, n, max(T(:, 2)));
                    T = obj.auxFlux.boundary.assembleJump(2*dim, cR*C, f, g, varargin{:});
                    obj.trJmpMat{2*dim} = core.linalg.sparseFromTriplets(T, n, max(T(:, 2)));
                end

                for i = 1:obj.space.nDims
                    for j = 1:obj.space.nDims
                        if coe(i, j) == 0, continue; end
                        switch obj.auxFlux.fluxType
                            case 'left'
                                J = obj.trJmpMat{2*j-1};
                            case 'right'
                                J = obj.trJmpMat{2*j};
                            case 'central'
                                J = (obj.trJmpMat{2*j-1} + obj.trJmpMat{2*j}) / 2;
                        end
                        A = A + coe(i, j) * J;
                    end
                end
            end
        end
    end
end