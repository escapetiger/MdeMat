classdef FluxAssembly < fem.assembly.Assembly
    % FLUXASSEMBLY Assembly for flux terms.

    properties
        bcType   % Boundary condition type: {'periodic', 'dirichlet'}
        fluxType % Flux type: {'left', 'right', 'central'}
        interior % Flux assembly on the interior elements
        boundary % Flux assembly on the boundary elements
    end

    methods
        function obj = FluxAssembly(context, bcType, fluxType)
            % FLUXASSEMBLY Constructor for FluxAssembly.
            %
            %   obj = FluxAssembly(context, bcType, fluxType)
            %   creates a flux assembly object with the specified boundary
            %   condition type and flux type.
            %
            % Inputs:
            %   context - Context object
            %   bcType - Boundary condition type: {'periodic', 'dirichlet'}
            %   fluxType - Flux type: {'left', 'right', 'central'}
            %   fluxMode - Flux mode: {'upwind', 'auxiliary', 'primal'}
            %
            % Outputs:
            %   obj - Constructed FluxAssembly object

            obj@fem.assembly.Assembly(context);
            obj.setBcType(bcType);
            obj.setFluxType(fluxType);
            obj.interior = fem.assembly.InteriorFluxAssembly(context);
            switch lower(obj.bcType)
                case 'periodic'
                    cls = 'PeriodicBoundaryFluxAssembly';
                case 'dirichlet'
                    cls = 'DirichletBoundaryFluxAssembly';
            end
            obj.boundary = fem.assembly.(cls)(context);
        end


        function obj = setFluxType(obj, fluxType)
            % SETFLUXTYPE Set the flux type.
            %
            %   obj = setFluxType(obj, fluxType) set the flux type.
            %
            % Inputs:
            %   obj - The FluxAssembly object
            %   fluxType - Flux type: {'left', 'right', 'central'}
            %
            % Outputs:
            %   obj - The FluxAssembly object

            core.except.assert(isempty(fluxType) || ...
                ismember(fluxType, {'left', 'right', 'central'}), ...
                'InvalidInput', 'Flux type is not supported.');

            obj.fluxType = fluxType;
        end

        function obj = setFluxMode(obj, fluxMode)
            % SETFLUXMODE Set the flux mode.
            %
            %   obj = setFluxMode(obj, fluxMode) set the flux mode.
            %
            % Inputs:
            %   obj - The FluxAssembly object
            %   fluxMode - Flux mode: {'upwind', 'auxiliary', 'primal'}
            %
            % Outputs:
            %   obj - The FluxAssembly object

            core.except.assert(strcmpi(obj.bcType, 'dirichlet'), ...
                'InvalidOperation', 'Flux mode requires Dirichlet BC.');

            obj.boundary.setFluxMode(fluxMode);
        end

        function T = assembleFluxPartial(obj, dim, coe)
            % ASSEMBLEFLUXPARTIAL Assemble flux terms for partial operator.
            %
            %   T = assembleFluxPartial(obj, dim, coe) creates triplets for
            %   the flux contribution of partial operators along the
            %   dimension @a dim with the coefficient @a coe.
            %
            % Inputs:
            %   obj - The FluxAssembly object
            %   dim - Dimension index
            %   coe - Coefficient vector
            %
            % Outputs:
            %   T - Triplet matrix [row, col, val] for flux terms
            switch lower(obj.fluxType)
                case 'left'
                    T = obj.assembleLeftFlux(dim, coe);
                case 'right'
                    T = obj.assembleRightFlux(dim, coe);
                case 'central'
                    T = obj.assembleCentralFlux(dim, coe);
            end
        end
    
        function T = assembleImplicitBoundaryJump(obj, dim, coe)
            % ASSEMBLEFLUXPARTIAL Assemble implicit boundary jump term.
            %
            %   T = assembleImplicitBoundaryJump(obj, dim, coe) creates
            %   triplets for the jump contribution of boundary flux along
            %   the dimension @a dim with the coefficient @a coe.
            %
            % Inputs:
            %   obj - The FluxAssembly object
            %   dim - Dimension index
            %   coe - Coefficient vector
            %
            % Outputs:
            %   T - Triplet matrix [row, col, val] for jump terms

            switch lower(obj.fluxType)
                case 'left'
                    T = obj.boundary.assembleImplicitLeftJump(dim, coe);
                case 'right'
                    T = obj.boundary.assembleImplicitRightJump(dim, coe);
                case 'central'
                    T1 = obj.boundary.assembleImplicitLeftJump(dim, coe/2);
                    T2 = obj.boundary.assembleImplicitRightJump(dim, coe/2);
                    T = [T1; T2];
            end
        end
    end

    methods (Access = protected)
        function obj = setBcType(obj, bcType)
            % SETBCTYPE Set the boundary condition type.
            %
            %   obj = setBcType(obj, bcType) set the boundary condition
            %   type.
            %
            % Inputs:
            %   obj - The FluxAssembly object
            %   bcType - Boundary condition type: {'periodic', 'dirichlet'}
            %
            % Outputs:
            %   obj - The FluxAssembly object

            core.except.assert(isempty(bcType) || ...
                ismember(bcType, {'periodic', 'dirichlet'}), ...
                'InvalidInput', 'Boundary condition type is not supported.');

            obj.bcType = bcType;
        end

        function T = assembleLeftFlux(obj, dim, coe)
            % ASSEMBLELEFTFLUX Assemble left-biased flux terms.
            %
            %   T = assembleLeftFlux(obj, dim, coe) combines interior and
            %   boundary contributions for left-biased flux terms along the
            %   specified dimension.
            %
            % Inputs:
            %   obj - The FluxAssembly object
            %   dim - Dimension index
            %   coe - Coefficient vector
            %
            % Outputs:
            %   T - Combined triplet matrix for left flux terms

            T1 = obj.interior.assembleInnerLeftFlux(dim, coe);
            T2 = obj.interior.assembleOuterLeftFlux(dim, coe);
            T3 = obj.boundary.assembleInnerLeftFlux(dim, coe);
            T4 = obj.boundary.assembleOuterLeftFlux(dim, coe);
            T = [T1; T2; T3; T4];
        end

        function T = assembleRightFlux(obj, dim, coe)
            % ASSEMBLYIGHTFLUX Assemble right-biased flux terms.
            %
            %   T = assembleRightFlux(obj, dim, coe) combines interior and
            %   boundary contributions for right-biased flux terms along
            %   the specified dimension.
            %
            % Inputs:
            %   obj - The FluxAssembly object
            %   dim - Dimension index
            %   coe - Coefficient vector
            %
            % Outputs:
            %   T - Combined triplet matrix for right flux terms

            T1 = obj.interior.assembleInnerRightFlux(dim, coe);
            T2 = obj.interior.assembleOuterRightFlux(dim, coe);
            T3 = obj.boundary.assembleInnerRightFlux(dim, coe);
            T4 = obj.boundary.assembleOuterRightFlux(dim, coe);
            T = [T1; T2; T3; T4];
        end

        function T = assembleCentralFlux(obj, dim, coe)
            % ASSEMBLECENTRALFLUX Assemble central flux terms.
            %
            %   T = assembleCentralFlux(obj, dim, coe) combines left and
            %   right flux contributions with equal weighting to create
            %   central difference-type operators.
            %
            % Inputs:
            %   obj - The FluxAssembly object
            %   dim - Dimension index
            %   coe - Coefficient vector
            %
            % Outputs:
            %   T - Combined triplet matrix for central flux terms

            T1 = obj.assembleLeftFlux(dim, coe/2);
            T2 = obj.assembleRightFlux(dim, coe/2);
            T = [T1; T2];
        end
    end
end
