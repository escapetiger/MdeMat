classdef HermiteState < physics.State
    % HERMITESTATE 
    
    properties
        nModes
    end
    
    methods
        function obj = HermiteState(xDisc, nModes)
            obj@physics.State(xDisc);
            obj.dofs = struct( ...
                'D', [], ... %< Hermite coefficients
                'omega', [] ... %< Modified potential
                );
            obj.coefs = struct( ...
                'EInf', [], ... %< Electrical field steady state
                'rhoInf', [] ... %< Density steady state
                );
            obj.nModes = nModes;
        end
    end
end

