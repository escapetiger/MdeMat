classdef GridAssembly < fem.assembly.Assembly
    % GRIDASSEMBLY Grid-based assembly for structured meshes.

    methods
        function obj = GridAssembly(fe, mesh, op)
            core.except.assert(isa(mesh, 'approx.mesh.Grid'), ...
                'InvalidInput', ...
                'Finite element space must be based on grid.')

            obj@fem.assembly.Assembly(fe, mesh, op);
        end
    end
end