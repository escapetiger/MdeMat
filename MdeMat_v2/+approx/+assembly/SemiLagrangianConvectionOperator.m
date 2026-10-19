classdef SemiLagrangianConvectionOperator < approx.assembly.Assembly
    % SEMILAGRANGIANCONVECTIONOPERATOR Assembly for semi-Lagrangian
    % convection operator.
    %
    %   SemiLagrangianConvectionOperator assembles linear transport
    %   operators for semi-Lagrangian methods. It handles characteristic
    %   tracing and flux reconstruction across element boundaries for
    %   advection-dominated problems, supporting both periodic and
    %   Dirichlet boundary conditions.
    %
    %   The operator uses backward characteristic tracing to compute the
    %   contribution from departure points to arrival points, enabling
    %   stable time stepping for convection problems without CFL
    %   restrictions. The assembly process coordinates volume integrals for
    %   interior contributions and boundary integrals for flux conditions.
    %
    % Examples:
    %   % Create semi-Lagrangian convection operator with periodic BCs
    %   space = approx.space.MeshSpace(element, mesh);
    %   slOperator = approx.element.L2SemiLagrangianElementOperator(element, clipper);
    %   convOp = SemiLagrangianConvectionOperator(space, slOperator, 'periodic');
    %   
    %   % Assemble for pure advection (periodic case)
    %   A = convOp.linear();
    %   
    %   % Assemble with Dirichlet boundary conditions
    %   convOpDir = SemiLagrangianConvectionOperator(space, slOperator, 'dirichlet');
    %   tStart = 0.0; tEnd = 0.1;
    %   bcFunc = @(x, t) sin(pi * x);
    %   [b, A] = convOpDir.linear(tStart, tEnd, bcFunc);
    %   
    %   % Use in time stepping
    %   U_new = A * U_old + b;
    %
    % Notes:
    %   Element operator must be a semi-Lagrangian type. Boundary condition
    %   type determines the assembly strategy: periodic uses volume
    %   integrals only, while Dirichlet includes boundary condition
    %   contributions.
    %
    % See also:
    %   approx.assembly.Assembly, 
    %   approx.element.SemiLagrangianElementOperator,
    %   approx.assembly.SemiLagrangianPeriodicAssembly,
    %   approx.assembly.SemiLagrangianDirichletAssembly

    properties
        dynamic % Dynamic assembly object
    end

    methods
        function obj = SemiLagrangianConvectionOperator(space, operator, bcType)
            % SEMILAGRANGIANCONVECTIONOPERATOR Constructor for
            % SemiLagrangianConvectionOperator.
            %
            %   obj = SemiLagrangianConvectionOperator(space, operator,
            %   bcType) creates a semi-Lagrangian convection operator with
            %   the specified function space, element operator, and
            %   boundary condition type. The constructor validates inputs
            %   and configures the appropriate assembly strategy.
            %
            % Inputs:
            %   space - Function space object (MeshSpace or similar)
            %   operator - Semi-Lagrangian element operator
            %   bcType - Boundary condition type: {'periodic', 'dirichlet'}
            %
            % Outputs:
            %   obj - Constructed SemiLagrangianConvectionOperator object

            core.except.assert( ...
                isa(operator, 'approx.element.SemiLagrangianElementOperator'), ...
                'SemiLagrangianConvectionOperator:InvalidOperator', ...
                'Element operator must be a Semi-Lagrangian element operator.');
                
            obj@approx.assembly.Assembly(space, operator);
            
            if strcmpi(bcType, 'periodic')
                obj.dynamic = approx.assembly.SemiLagrangianPeriodicAssembly(space, operator);
            else
                obj.dynamic = approx.assembly.SemiLagrangianDirichletAssembly(space, operator);
            end
        end

        function [b, A] = linear(obj, varargin)
            % LINEAR Assemble linear transport operator components.
            %
            %   b = linear(obj) assembles boundary condition vector for the
            %   linear transport operator with current boundary conditions.
            %   For periodic boundary conditions, returns empty vector.
            %
            %   [b, A] = linear(obj) assembles both the boundary condition
            %   vector and the linear transport operator matrix. The matrix
            %   represents the semi-Lagrangian discretization of the advection
            %   operator using characteristic tracing.
            %
            %   [b, A] = linear(obj, tBegin, tEnd, f) assembles linear
            %   transport components with specified time interval and
            %   boundary condition function. Time interval [tBegin, tEnd]
            %   defines the characteristic tracing period.
            %
            % Inputs:
            %   obj - The SemiLagrangianConvectionOperator object
            %   varargin - Variable input arguments for Dirichlet conditions:
            %<    tBegin - Beginning time for characteristic tracking (scalar)
            %<    tEnd - End time for characteristic tracking (scalar)
            %<    f - Boundary condition function handle f(x, t)
            %
            % Outputs:
            %   b - Boundary contribution vector (sparse column vector)
            %   A - Interior contribution matrix (sparse matrix)

            n = obj.space.nGlobalDofs;

            %< Interior contribution via characteristic tracing
            if nargout >= 2
                T = obj.dynamic.assembleVolume();
                A = core.linalg.sparseFromTriplets(T, n, n);
            end

            %< Boundary contribution for Dirichlet conditions
            if nargout >= 1
                b = [];
                if strcmpi(obj.dynamic.BC_TYPE, 'dirichlet')
                    [tBegin, tEnd, f] = varargin{:};
                    if abs(tBegin-tEnd) < 1e-8, return; end
                    T = obj.dynamic.assembleBoundaryCondition(tBegin, tEnd, f);
                    b = core.linalg.sparseFromTriplets(T, n, 1);
                end
            end
        end
    end
end