classdef State < handle
    % STATE Solution state for advection problems.
    %
    %   The State class encapsulates the complete solution state for
    %   advection problems, managing both spatial discretization
    %   information and degrees of freedom data. It provides a unified
    %   container for all components needed to represent and evolve
    %   the numerical solution.
    %
    %   The class uses handle semantics to enable efficient state
    %   management during time integration without unnecessary copying
    %   of large data structures.
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
    %
    % See Also:
    %   physics.Scheme, physics.advection.SldgScheme,
    %   physics.advection.UwdgScheme

    properties
        disc % Discretization
        dofs % Degrees of freedom
        coefs % Coefficients
    end

    methods
        function obj = State(disc)
            % STATE Constructor for State.
            %
            %   obj = State(disc) creates a solution state object with
            %   the specified discretization and initializes
            %   empty degrees of freedom structure.
            %
            % Inputs:
            %   disc - Discretization structure 
            %
            % Outputs:
            %   obj - Constructed State object

            obj.disc = disc;
            obj.dofs = struct('U', []);
            obj.coefs = struct();
        end
    end
end