classdef State < handle
    % STATE Base class for all solution states.
    %
    %   The State class encapsulates the complete solution state, managing
    %   discretization, degrees of freedom and coefficients. It provides a
    %   unified container for all components needed to represent and evolve
    %   the numerical solution.
    %
    % Examples:
    %   % Create state with spatial discretization
    %   state = State(spatialDiscretization);
    %   
    %   % Access discretization and solution data
    %   mesh = state.disc.x.mesh;
    %   solution = state.dofs.U;
    %   
    %   % Update solution data
    %   state.dofs.U = newSolutionVector;
    %
    % Notes:
    %   The handle-based design allows for efficient passing of state
    %   objects between methods without copying large arrays.

    properties
        xDisc % Spatial discretization
        dofs % Degrees of freedom
        coefs % Coefficients
    end

    methods
        function obj = State(xDisc)
            % STATE Constructor for State.
            %
            %   obj = State(xDisc) creates a solution state object with the
            %   specified spatial discretization and initializes empty
            %   degrees of freedom structure.
            %
            % Inputs:
            %   xDisc - FiniteElementDiscretization structure 
            %
            % Outputs:
            %   obj - Constructed State object

            obj.xDisc = xDisc;
            obj.dofs = struct('U', []);
            obj.coefs = struct();
        end
    end
end