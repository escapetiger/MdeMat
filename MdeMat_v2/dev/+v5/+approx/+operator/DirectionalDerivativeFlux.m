classdef DirectionalDerivativeFlux < approx.operator.MeshOperator
    properties
        VALID_FLUX_TYPES = {'upwind', 'central', 'outward', 'inward'}
        VALID_PERIODIC_FLUX_MODES = {'periodic'}
        VALID_DIRICHLET_FLUX_MODES = {'extrapolate', 'truncate', 'inflow'}
        
    end
    
    properties
        bcType % Boundary condition type
        fluxType % Flux type
        fluxMode % Flux mode
    end

    methods

    end
end

