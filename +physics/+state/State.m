classdef State < handle
    % STATE Base class for all solution states.
    %
    %   State encapsulates the fundamental solution state components
    %   including degrees of freedom and coefficients. It provides
    %   a unified container and interface for state management operations
    %   including reporting, saving, and loading.
    %
    %   The handle-based design allows for efficient passing of state
    %   objects between methods without copying large arrays.
    %
    % See also:
    %   physics.state.KineticState, physics.state.SpatialState
    
    properties
        Dofs    % Degrees of freedom structure
        Coefs   % Coefficients structure
        History % Time series data structure
    end
    
    methods
        function obj = State()
            % STATE Construct an instance of State.
            %
            %   obj = State() creates a base solution state object with
            %   empty degrees of freedom and coefficients structures.
            
            obj.Dofs = struct();
            obj.Coefs = struct();
            obj.History = struct();
        end
        
        function obj = setDof(obj, name, value)
            % SETDOF Set degree of freedom in state.
            %
            %   obj = setDof(obj, name, value) sets a named degree of
            %   freedom with the specified @a value in the state.
            
            arguments
                obj physics.state.State
                name {mustBeTextScalar}
                value {mustBeNumeric} = []
            end
            
            obj.Dofs.(char(name)) = value;
        end
        
        function obj = setCoefficient(obj, name, value)
            % SETCOEFFICIENT Set coefficient in state.
            %
            %   setCoefficient(obj, name, value) sets a named coefficient
            %   with the specified @a value in the state.
            
            arguments
                obj physics.state.State
                name {mustBeTextScalar}
                value {mustBeNumeric} = []
            end
            
            obj.Coefs.(char(name)) = value;
        end

        function obj = setHistory(obj, name, value)
            % SETHISTORY Set time series data in state history.
            %
            %   obj = setHistory(obj, name, value) sets a named time
            %   series @a value in the state's history.

            arguments
                obj physics.state.State
                name {mustBeTextScalar}
                value {mustBeNumeric} = []
            end

            obj.History.(char(name)) = value;
        end
        
        function obj = report(obj)
            % REPORT Generate state information report.
            %
            %   obj = report(obj) displays the names of degrees of freedom
            %   and coefficients stored in the state.

            arguments
                obj physics.state.State
            end
            
            fprintf('[R] State Type: %s\n', class(obj));
            
            dofNames = fieldnames(obj.Dofs);
            if isempty(dofNames)
                fprintf('[R] DOF Names: (none)\n');
            else
                fprintf('[R] DOF Names: %s\n', strjoin(dofNames, ', '));
            end
            
            coefNames = fieldnames(obj.Coefs);
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
            %   to a MAT file with the specified @a filename. Automatically
            %   creates the directory if it does not exist.
            
            arguments
                obj physics.state.State
                filename {mustBeTextScalar}
            end
            
            [dirPath, ~, ~] = fileparts(filename);
            if ~isempty(dirPath) && ~exist(dirPath, 'dir')
                mkdir(dirPath);
                fprintf('[M] Created directory: %s\n', dirPath);
            end
            
            state = struct();
            state.dofs = obj.Dofs;
            state.coefs = obj.Coefs;
            state.history = obj.History;
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
            
            arguments
                filename {mustBeTextScalar}
            end
            
            core.except.assert(exist(filename, 'file') == 2, ...
                'InvalidInput', 'File %s does not exist', filename);
            
            state = load(filename).state;
            
            if isfield(state, 'className')
                obj = feval(state.className);
            else
                obj = physics.state.State();
            end
            
            obj.Dofs = state.dofs;
            obj.Coefs = state.coefs;
            obj.History = state.history;

            fprintf('[M] State loaded from: %s\n', filename);
        end
    end
end