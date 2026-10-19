classdef AuxiliaryDivergenceOperator < approx.assembly.Assembly
    % AUXILIARYDIVERGENCEOPERATOR Assembly for auxiliary divergence
    % operator.

    properties
        volume % Volume assembly for interior terms (VolumeAssembly)
        flux % Auxiliary flux assembly (FluxAssembly)
        auxMat % Auxiliary operator matrices (cell array)
        jmpMat % Implicit jump operator matrices (cell array)
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
        function obj = AuxiliaryDivergenceOperator(space, operator, bcType, fluxType, jumpType)
            % DIFFUSIONOPERATOR Constructor for AuxiliaryDivergenceOperator.
            %
            %   obj = AuxiliaryDivergenceOperator(space, operator, bcType,
            %   fluxType, jumpType) creates an auxiliary divergence
            %   operator with the specified components for Local
            %   Discontinuous Galerkin discretization.
            %
            % Inputs:
            %   space - MeshSpace object
            %   operator - ElementOperator object
            %   bcType - Boundary condition type: {'periodic', 'dirichlet'}
            %   fluxType - Flux type: {'left', 'right', 'central'}
            %   jumpType - Jump treatment: {'implicit', 'explicit'}
            %
            % Outputs:
            %   obj - Constructed AuxiliaryDivergenceOperator object

            obj@approx.assembly.Assembly(space, operator);
            obj.volume = approx.assembly.VolumeAssembly(space, operator);
            obj.flux = approx.assembly.FluxAssembly(space, operator, bcType, fluxType);
            if strcmpi(bcType, 'dirichlet')
                obj.flux.setFluxMode('auxiliary');
            end
            obj.auxMat = cell(1, obj.space.nDims);
            obj.jmpMat = cell(1, obj.space.nDims);
            obj.trMat = cell(1, obj.space.element.nFluxes);
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

            TF = strcmpi(obj.flux.bcType, 'dirichlet') && strcmpi(obj.jumpType, 'implicit');
        end

        function TF = get.isExplicitDirichletBoundaryJump(obj)
            % GET.ISEXPLICITDIRICHLETBOUNDARYJUMP Check for explicit
            % Dirichlet jump treatment.

            TF = strcmpi(obj.flux.bcType, 'dirichlet') && strcmpi(obj.jumpType, 'explicit');
        end

        function A = linear(obj, coe)
            % LINEAR Assemble auxiliary linear divergence matrix.
            %
            %   A = linear(obj, coe) assembles the auxiliary divergence of
            %   linear field.
            %
            % Inputs:
            %   obj - The AuxiliaryDivergenceOperator object
            %   coe - Constant coefficients (nDims x 1 vector)
            %
            % Outputs:
            %   A - Assembled matrix (sparse)

            n = obj.space.nGlobalDofs;
            h = obj.space.mesh.spacings;

            for dim = 1:obj.space.nDims
                if isempty(obj.auxMat{dim})
                    C = obj.volume.scaleConstant(dim, 1);
                    T1 = obj.volume.assembleVolumePartial(dim, C);
                    T2 = obj.flux.assembleFluxPartial(dim, C);
                    T = [T1; T2];
                    obj.auxMat{dim} = core.linalg.sparseFromTriplets(T, n, n);
                end
            end

            A = sparse(n, n);
            for dim = 1:obj.space.nDims
                if coe(dim) == 0, continue; end
                A = A + coe(dim) * obj.auxMat{dim};
            end
        end

        function A = linearJumpExtraction(obj, coe)
            % LINEARJUMPEXTRACTION Assemble the extraction part of jump
            % term for auxiliary linear divergence operator.
            %
            %   A = linearJump(obj, coe) assembles the extraction part of
            %   implicit jump term for auxiliary divergence of linear
            %   field.
            %
            % Inputs:
            %   obj - The AuxiliaryDivergenceOperator object
            %   coe - Constant coefficients (nDims x 1 vector)
            %
            % Outputs:
            %   A - Assembled matrix (sparse)

            n = obj.space.nGlobalDofs;
            h = obj.space.mesh.spacings;

            if obj.isImplicitDirichletBoundaryJump
                for dim = 1:obj.space.nDims
                    if isempty(obj.jmpMat{dim})
                        C = obj.volume.scaleConstant(dim, 1);
                        switch lower(obj.flux.fluxType)
                            case 'left'
                                T = obj.flux.boundary.assembleImplicitLeftJump(dim, C);
                            case 'right'
                                T = obj.flux.boundary.assembleImplicitRightJump(dim, C);
                            case 'central'
                                T1 = obj.flux.boundary.assembleImplicitLeftJump(dim, C/2);
                                T2 = obj.flux.boundary.assembleImplicitRightJump(dim, C/2);
                                T = [T1; T2];
                        end
                        obj.jmpMat{dim} = core.linalg.sparseFromTriplets(T, n, n);
                    end
                end

                A = sparse(n, n);
                for dim = 1:obj.space.nDims
                    if coe(dim) == 0, continue; end
                    A = A + coe(dim) * obj.jmpMat{dim};
                end
            end
        end
    end
end