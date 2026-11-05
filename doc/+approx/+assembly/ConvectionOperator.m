classdef ConvectionOperator < approx.assembly.Assembly
    % CONVECTIONOPERATOR Assembly for convection operators in DG methods.
    %
    %   ConvectionOperator provides high-level assembly of convection
    %   (advection) operators using discontinuous Galerkin discretizations.
    %   It combines volume and flux assembly components to create complete
    %   convection operators with appropriate upwind stabilization.
    %
    %   The operator automatically selects upwind flux directions based
    %   on the velocity field, ensuring stability for convection-dominated
    %   problems. It supports both periodic and Dirichlet boundary
    %   conditions with appropriate flux treatments.
    %
    %   The convection operator discretizes terms of the form:
    %
    %   \f[
    %     \nabla \cdot (c \mathbf{v} u)
    %   \f] 
    % 
    %   where \f$c\f$ is a coefficient, \f$v\f$ is the velocity field,
    %   and \f$u\f$ is the transported quantity.
    %
    % Examples:
    %   % Create convection operator with periodic boundaries
    %   space = approx.space.MeshSpace(mesh, element);
    %   operator = approx.element.L2ElementOperator(element);
    %   convOp = ConvectionOperator(space, operator, 'periodic');
    %
    %   % Assemble convection matrix for constant velocity
    %   coefficients = [1.0; 0.5];  % [c_x; c_y]
    %   velocity = [2.0; 1.0];      % [v_x; v_y]
    %   A = convOp.linear(coefficients, velocity);
    %
    %   % Assemble boundary condition terms
    %   boundaryFunction = @(x, t) sin(pi*x(1,:)) .* cos(pi*x(2,:));
    %   B = convOp.linearBc(coefficients, velocity, boundaryFunction);
    %
    % See also:
    %   approx.assembly.Assembly, approx.assembly.VolumeAssembly,
    %   approx.assembly.FluxAssembly

    properties
        volume % Volume assembly for interior terms (VolumeAssembly)
        flux   % Flux assembly for interface terms (FluxAssembly)
    end

    methods
        function obj = ConvectionOperator(space, operator, bcType)
            % CONVECTIONOPERATOR Constructor for ConvectionOperator.
            %
            %   obj = ConvectionOperator(space, operator, bcType) creates a
            %   convection operator with the specified mesh space, element
            %   operator, and boundary condition type. For Dirichlet
            %   boundaries, upwind flux is automatically selected.
            %
            % Inputs:
            %   space - MeshSpace object
            %   operator - ElementOperator object
            %   bcType - Boundary condition type: {'periodic', 'dirichlet'}
            %
            % Outputs:
            %   obj - Constructed ConvectionOperator object

            obj@approx.assembly.Assembly(space, operator);
            obj.volume = approx.assembly.VolumeAssembly(space, operator);
            obj.flux = approx.assembly.FluxAssembly(space, operator, bcType, []);
            if strcmpi(bcType, 'dirichlet')
                obj.flux.setFluxMode('upwind');
            end
        end

        function A = linear(obj, coe, vel)
            % LINEAR Assemble linear convection operator matrix.
            %
            %   A = linear(obj, coe) assembles the convection operator
            %   matrix with coefficients coe, using coe as both
            %   coefficients and velocity field.
            %
            %   A = linear(obj, coe, vel) assembles the convection operator
            %   matrix with coefficients coe and velocity field vel. The
            %   flux direction (left/right bias) is automatically selected
            %   based on the velocity sign in each dimension.
            %
            % Inputs:
            %   obj - The ConvectionOperator object
            %   coe - Coefficient vector for each dimension (nDims x 1)
            %   vel - Velocity vector for each dimension (nDims x 1, optional)
            %
            % Outputs:
            %   A - Assembled convection operator matrix (sparse)

            if nargin < 3 || isempty(vel), vel = coe; end

            n = obj.space.nGlobalDofs;
            A = sparse(n, n);
            for dim = 1:obj.space.nDims
                if coe(dim) == 0, continue; end
                if vel(dim) > 0
                    obj.flux.setFluxType('left');
                else
                    obj.flux.setFluxType('right');
                end
                C = obj.volume.scaleConstant(dim, coe(dim));
                T1 = obj.volume.assembleVolumePartial(dim, C);
                T2 = obj.flux.assembleFluxPartial(dim, C);
                T = [T1; T2];
                B = core.linalg.sparseFromTriplets(T, n, n);
                A = A + B;
            end
        end

        function A = linearBc(obj, coe, vel, f, varargin)
            % LINEARBC Assemble boundary condition terms for linear
            % convection.
            %
            %   A = linearBc(obj, coe, vel, f) assembles the boundary
            %   condition contribution matrix for linear convection with
            %   Dirichlet boundary conditions. The boundary flux direction
            %   depends on the velocity sign (inflow boundaries only).
            %
            % Inputs:
            %   obj - The ConvectionOperator object
            %   coe - Coefficient vector for each dimension (nDims x 1)
            %   vel - Velocity vector for each dimension (nDims x 1)
            %   f - Boundary condition function handle
            %   varargin - Input arguments passed to boundary condition
            %
            % Outputs:
            %   A - Boundary condition contribution matrix (sparse)

            n = obj.space.nGlobalDofs;
            A = cell(1, obj.space.nDims);
            for dim = 1:obj.space.nDims
                if coe(dim) == 0, continue; end
                C = obj.volume.scaleConstant(dim, coe(dim));
                if vel(dim) > 0
                    T = obj.flux.boundary.assembleTrace(2*dim-1, -C, f, varargin{:});
                else
                    T = obj.flux.boundary.assembleTrace(2*dim, C, f, varargin{:});
                end
                A{dim} = core.linalg.sparseFromTriplets(T, n, max(T(:, 2)));
            end

            for dim = obj.space.nDims:-1:2
                A{1} = A{1} + A{dim};
                A{dim} = [];
            end
            A = A{1};
        end
    end
end