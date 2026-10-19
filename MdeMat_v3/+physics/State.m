classdef State < handle
    % STATE Base class for all solution states.
    %
    %   The State class encapsulates the fundamental solution state
    %   components including degrees of freedom and coefficients. It
    %   provides a unified container and interface for state management
    %   operations including reporting, saving, and loading.
    %
    % Examples:
    %   % Create basic state
    %   state = physics.State();
    %   
    %   % Save and load operations
    %   state.save('solution.mat');
    %   newState = state.load('solution.mat');
    %
    % Notes:
    %   The handle-based design allows for efficient passing of state
    %   objects between methods without copying large arrays.

    properties
        dofs %< Degrees of freedom structure
        coefs %< Coefficients structure
    end

    methods
        function obj = State()
            % STATE Constructor for State.
            %
            %   obj = State() creates a base solution state object with
            %   empty degrees of freedom and coefficients structures.
            %
            % Outputs:
            %   obj - Constructed State object

            obj.dofs = struct();
            obj.coefs = struct();
        end

        function obj = addDof(obj, name, value)
            % ADDDOF Add degree of freedom to state.
            %
            %   obj = addDof(obj, name, value) adds a named degree of
            %   freedom with the specified value to the state. Value can be
            %   empty.
            %
            % Inputs:
            %   obj - The State object
            %   name - DOF name (char or string)
            %   value - DOF value (any type, can be [])
            %
            % Outputs:
            %   obj - The State object

            obj.dofs.(char(name)) = value;
        end

        function addCoefficient(obj, name, value)
            % ADDCOEFFICIENT Add coefficient to state.
            %
            %   obj = addCoefficient(obj, name, value) adds a named
            %   coefficient with the specified value to the state. Value
            %   can be empty.
            %
            % Inputs:
            %   obj - The State object
            %   name - Coefficient name (char or string)
            %   value - Coefficient value (any type, can be [])
            %
            % Outputs:
            %   obj - The State object

            obj.coefs.(char(name)) = value;
        end

        function obj = report(obj)
            % REPORT Generate state information report.
            %
            %   obj = report(obj) displays the names of degrees of freedom
            %   and coefficients stored in the state.
            %
            % Inputs:
            %   obj - The State object
            %
            % Outputs:
            %   obj - The State object

            fprintf('[R] State Type: %s\n', class(obj));
            
            %< Report DOF names
            dofNames = fieldnames(obj.dofs);
            if isempty(dofNames)
                fprintf('[R] DOF Names: (none)\n');
            else
                fprintf('[R] DOF Names: %s\n', strjoin(dofNames, ', '));
            end
            
            %< Report coefficient names
            coefNames = fieldnames(obj.coefs);
            if isempty(coefNames)
                fprintf('[R] Coefficient Names: (none)\n');
            else
                fprintf('[R] Coefficient Names: %s\n', strjoin(coefNames, ', '));
            end
        end

        function obj = save(obj, filename)
            % SAVE Save state to file.
            %
            %   obj = save(obj, filename) saves the complete state object
            %   to a MAT file with the specified filename.
            %
            % Inputs:
            %   obj - The State object
            %   filename - Output filename (char or string)
            %
            % Outputs:
            %   obj - The State object

            %< Input validation (rule 235)
            core.except.assert(ischar(filename) || isstring(filename), ...
                'InvalidInput', 'filename must be char or string');

            %< Save state data
            state.dofs = obj.dofs;
            state.coefs = obj.coefs;
            state.className = class(obj);
            
            save(filename, 'state', '-v7.3');
            
            fprintf('[M] State saved to: %s\n', filename);
        end
    end

    methods (Static)
        function obj = load(filename)
            % LOAD Load state from file.
            %
            %   obj = State.load(filename) loads a state object from the
            %   specified MAT file and returns a new instance of the
            %   appropriate class.
            %
            % Inputs:
            %   filename - Input filename (char or string)
            %
            % Outputs:
            %   obj - New State object with loaded data

            %< Input validation (rule 235)
            core.except.assert(ischar(filename) || isstring(filename), ...
                'InvalidInput', 'filename must be char or string');
            core.except.assert(exist(filename, 'file') == 2, ...
                'InvalidInput', 'File %s does not exist', filename);

            %< Load state data
            state = load(filename).state;
            
            %< Create new object of stored class type
            if isfield(state, 'className')
                obj = feval(state.className);
            else
                obj = State();
            end
            
            %< Restore state data
            obj.dofs = state.dofs;
            obj.coefs = state.coefs;
            
            fprintf('[M] State loaded from: %s\n', filename);
        end
    end
end