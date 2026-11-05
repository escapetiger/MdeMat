classdef FluxAssembly < approx.assembly.Assembly
    % FLUXASSEMBLY Assembly for flux terms in DG methods.
    %
    %   FluxAssembly handles the assembly of flux terms that arise in
    %   discontinuous Galerkin (DG) formulations. These flux terms
    %   provide inter-element coupling in DG methods, where the solution
    %   is allowed to be discontinuous across element boundaries.
    %
    %   The class supports different flux types (left-biased, right-biased,
    %   central) and boundary condition types (periodic, Dirichlet).
    %   It manages both interior flux assembly (between neighboring
    %   elements) and boundary flux assembly (at domain boundaries).
    %
    %   Flux terms are essential for stability and accuracy in DG methods,
    %   as they control the flow of information across element interfaces
    %   and enforce appropriate boundary conditions. The choice of flux
    %   type affects the stability properties and numerical diffusion
    %   characteristics of the resulting scheme.
    %
    % See also:
    %   approx.assembly.Assembly, approx.assembly.VolumeAssembly,
    %   approx.assembly.InteriorFluxAssembly

    properties
        bcType   % Boundary condition type: {'periodic', 'dirichlet'}
        fluxType % Flux type: {'left', 'right', 'central'}
        interior % Flux assembly on the interior elements (InteriorFluxAssembly)
        boundary % Flux assembly on the boundary elements (BoundaryFluxAssembly)
    end

    methods
        function obj = FluxAssembly(space, operator, bcType, fluxType)
            % FLUXASSEMBLY Constructor for FluxAssembly.
            %
            %   obj = FluxAssembly(space, operator, bcType, fluxType)
            %   creates a flux assembly object with the specified boundary
            %   condition type and flux type. The constructor automatically
            %   creates appropriate interior and boundary flux assembly
            %   objects based on the specified boundary condition type.
            %
            % Inputs:
            %   space - MeshSpace object
            %   operator - ElementOperator object
            %   bcType - Boundary condition type: {'periodic', 'dirichlet'}
            %   fluxType - Flux type: {'left', 'right', 'central'}
            %
            % Outputs:
            %   obj - Constructed FluxAssembly object

            obj@approx.assembly.Assembly(space, operator);
            obj.setBcType(bcType);
            obj.setFluxType(fluxType);
            obj.interior = approx.assembly.InteriorFluxAssembly(space, operator);
            switch lower(obj.bcType)
                case 'periodic'
                    cls = 'PeriodicBoundaryFluxAssembly';
                case 'dirichlet'
                    cls = 'DirichletBoundaryFluxAssembly';
            end
            obj.boundary = approx.assembly.(cls)(space, operator);
        end

        function obj = setFluxType(obj, fluxType)
            % SETFLUXTYPE Set the flux type.
            %
            %   obj = setFluxType(obj, fluxType) sets the flux type for
            %   the assembly. The flux type determines how numerical
            %   fluxes are computed at element interfaces, affecting
            %   stability and accuracy properties.
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
            % SETFLUXMODE Set the flux mode for Dirichlet boundaries.
            %
            %   obj = setFluxMode(obj, fluxMode) sets the flux mode for
            %   Dirichlet boundary conditions. This method is only
            %   applicable when the boundary condition type is 'dirichlet'.
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
            %   dimension dim with the coefficient coe. The specific flux
            %   implementation depends on the flux type (left, right,
            %   central).
            %
            % Inputs:
            %   obj - The FluxAssembly object
            %   dim - Dimension index (positive integer)
            %   coe - Coefficient vector (nElements x 1)
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
            % ASSEMBLEIMPLICITBOUNDARYJUMP Assemble implicit boundary jump terms.
            %
            %   T = assembleImplicitBoundaryJump(obj, dim, coe) creates
            %   triplets for the jump contribution of boundary flux along
            %   the dimension dim with the coefficient coe. Jump terms
            %   arise in DG formulations to handle discontinuities at
            %   domain boundaries.
            %
            % Inputs:
            %   obj - The FluxAssembly object
            %   dim - Dimension index (positive integer)
            %   coe - Coefficient vector (nElements x 1)
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
            %   obj = setBcType(obj, bcType) sets the boundary condition
            %   type for the assembly. This determines how boundary
            %   fluxes are handled and which boundary assembly class
            %   is instantiated.
            %
            % Inputs:
            %   obj - The FluxAssembly object
            %   bcType - Boundary condition type: {'periodic', 'dirichlet'}
            %
            % Outputs:
            %   obj - The FluxAssembly object

            core.except.assert(isempty(bcType) || ismember(bcType, {'periodic', 'dirichlet'}), ...
                'InvalidInput', 'Boundary condition type is not supported.');

            obj.bcType = bcType;
        end

        function T = assembleLeftFlux(obj, dim, coe)
            % ASSEMBLELEFTFLUX Assemble left-biased flux terms.
            %
            %   T = assembleLeftFlux(obj, dim, coe) combines interior and
            %   boundary contributions for left-biased flux terms along the
            %   specified dimension. Left-biased fluxes provide upwind-like
            %   characteristics that can enhance stability for convection-
            %   dominated problems.
            %
            % Inputs:
            %   obj - The FluxAssembly object
            %   dim - Dimension index (positive integer)
            %   coe - Coefficient vector (nElements x 1)
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
            % ASSEMBLERIGHTFLUX Assemble right-biased flux terms.
            %
            %   T = assembleRightFlux(obj, dim, coe) combines interior and
            %   boundary contributions for right-biased flux terms along
            %   the specified dimension. Right-biased fluxes provide
            %   downwind-like characteristics.
            %
            % Inputs:
            %   obj - The FluxAssembly object
            %   dim - Dimension index (positive integer)
            %   coe - Coefficient vector (nElements x 1)
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
            %   central difference-type operators. Central fluxes provide
            %   symmetric treatment but may require additional
            %   stabilization for convection-dominated problems.
            %
            % Inputs:
            %   obj - The FluxAssembly object
            %   dim - Dimension index (positive integer)
            %   coe - Coefficient vector (nElements x 1)
            %
            % Outputs:
            %   T - Combined triplet matrix for central flux terms

            T1 = obj.assembleLeftFlux(dim, coe/2);
            T2 = obj.assembleRightFlux(dim, coe/2);
            T = [T1; T2];
        end
    end
end