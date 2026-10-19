classdef MixedHessian < fem.assembly.Operator    
    properties
        aux      % Auxiliary variable assembly structure
        prm      % Primal variable assembly structure
        jumpType % Jump treatment type (1=implicit, 2=explicit)
    end
    
    methods
        function obj = MixedHessian(fe, mesh, op, trace, extraction, auxFlux, prmFlux, jumpType)
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

            interior = fem.assembly.GridInteriorAssembly(fe, mesh, op);
            if isempty(trace)
                auxBoundary = fem.assembly.GridPeriodicBoundaryAssembly(fe, mesh, op);
            else
                auxBoundary = fem.assembly.AuxiliaryGridDirichletBoundaryAssembly(fe, mesh, op, trace, extraction);
            end
            if isempty(trace)
                prmBoundary = fem.assembly.GridPeriodicBoundaryAssembly(fe, mesh, op);
            else
                prmBoundary = fem.assembly.PrimalGridDirichletBoundaryAssembly(fe, mesh, op, trace);
            end

            obj@fem.assembly.Operator(fe, mesh, op);
            obj.aux.asm = fem.assembly.PartialGridAssembly(fe, mesh, op, interior, auxBoundary);
            obj.aux.mat = cell(1, obj.nDims);
            obj.prm.asm = fem.assembly.PrimalGridAssembly(fe, mesh, op, trace, extraction);
            obj.prm.mat = cell(1, obj.nDims);
            obj.jumpType = jumpType;
        end

    end
end

