classdef DiffusionOperator < fem.assembly.Assembly
    properties
        volume
        auxFlux
        prmFlux
        auxMat
        prmMat
        jmpMat
        jumpType % Jump treatment type (1=implicit, 2=explicit)
    end

    properties (Dependent)
        isImplicitDirichletBoundaryJump
        isExplicitDirichletBoundaryJump
    end

    methods
        function obj = DiffusionOperator(context, bcType, auxFluxType, prmFluxType, jumpType)
            obj.volume = fem.assembly.VolumeAssembly(context);
            obj.auxFlux = fem.assembly.FluxAssembly(context, bcType, auxFluxType);
            obj.prmFlux = fem.assembly.FluxAssembly(context, bcType, prmFluxType);
            if strcmpi(bcType, 'dirichlet')
                obj.auxFlux.setFluxMode('auxiliary');
                obj.prmFlux.setFluxMode('primal');
            end
            obj.auxMat = cell(1, obj.context.nDims);
            obj.prmMat = cell(1, obj.context.nDims);
            obj.jmpMat = cell(1, obj.context.nDims);
            obj.jumpType = jumpType;
        end

        function TF = get.isImplicitDirichletBoundaryJump(obj)
            TF = strcmpi(obj.auxFlux.bcType, 'dirichlet') && strcmpi(obj.jumpType, 'implicit');
        end

        function TF = get.isExplicitDirichletBoundaryJump(obj)
            TF = strcmpi(obj.auxFlux.bcType, 'dirichlet') && strcmpi(obj.jumpType, 'explicit');
        end

        function A = linear(obj, coe)
            n = obj.context.nGlobalDofs;

            for dim = 1:obj.context.nDims
                C = obj.volume.scaleConstant(dim, 1);
                T1 = obj.volume.assembleVolumePartial(dim, C);
                T2 = obj.auxFlux.assembleFluxPartial(dim, C);
                T = [T1; T2];
                obj.auxMat{dim} = core.linalg.sparseFromTriplets(T, n, n);

                T2 = obj.prmFlux.assembleFluxPartial(dim, C);
                T = [T1; T2];
                obj.prmMat{dim} = core.linalg.sparseFromTriplets(T, n, n);

                if obj.isImplicitDirichletBoundaryJump
                    T = obj.auxFlux.assembleImplicitBoundaryJump(dim, 1);
                    obj.jmpMat{dim} = core.linalg.sparseFromTriplets(T, n, n);
                end
            end

            A = sparse(n, n);
            for i = 1:obj.context.nDims
                for j = 1:obj.context.nDims
                    if coe(i, j) ~= 0
                        A = A + coe(i, j) * obj.auxMat{i} * obj.prmMat{j};
                        if obj.isImplicitDirichletBoundaryJump
                            A = A + coe(i, j) * obj.jmpMat{i};
                        end
                    end
                end
            end
        end

        function A = linearBc(obj, coe, f, g)
            n = obj.context.nGlobalDofs;
            A = cell(1, obj.context.nDims);
            for dim = 1:obj.context.nDims
                if coe(dim) == 0, continue; end
                C = obj.volume.assembleConstant(dim, coe(dim));
                T1 = obj.flux.boundary.assembleLeftTrace(dim, C, f);
                T2 = obj.flux.boundary.assembleRightTrace(dim, C, f);
                T = [T1; T2];
                A{dim} = core.linalg.sparseFromTriplets(T, n, max(T(:, 2)));
                A{dim} = obj.auxMat{dim} * A{dim};

                if obj.isImplicitDirichletBoundaryJump
                    switch obj.auxFlux.fluxType
                        case 'left'
                            T = T1;
                            B = core.linalg.sparseFromTriplets(T, n, max(T(:, 2)));
                        case 'right'
                            T = T2;
                            B = core.linalg.sparseFromTriplets(T, n, max(T(:, 2)));
                        case 'central'
                            T1(:, 3) = T1(:, 3) / 2;
                            T2(:, 3) = T2(:, 3) / 2;
                            T = [T1; T2];
                            B = core.linalg.sparseFromTriplets(T, n, max(T(:, 2)));
                    end
                    A{dim} = A{dim} + B;
                elseif obj.isExplicitDirichletBoundaryJump
                    switch obj.auxFlux.fluxType
                        case 'left'
                            T = obj.flux.boundary.assembleLeftJump(dim, C, f, g);
                            B = core.linalg.sparseFromTriplets(T, n, max(T(:, 2)));
                        case 'right'
                            T = obj.flux.boundary.assembleRightJump(dim, C, f, g);
                            B = core.linalg.sparseFromTriplets(T, n, max(T(:, 2)));
                        case 'central'
                            T1 = obj.flux.boundary.assembleLeftJump(dim, C / 2, f, g);
                            T2 = obj.flux.boundary.assembleRightJump(dim, C / 2, f, g);
                            T = [T1; T2];
                            B = core.linalg.sparseFromTriplets(T, n, max(T(:, 2)));
                    end
                    A{dim} = A{dim} + B;
                end
            end

            for dim = obj.context.nDims:-1:2
                A{1} = A{1} + A{dim};
                A{dim} = [];
            end
            A = A{1};
        end
    end
end
