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
    %   mdemat.common.fem.KineticState, mdemat.common.fem.SpatialState
    
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
                obj mdemat.data.variable.State
                name {mustBeTextScalar}
                value {mustBeNumeric} = []
            end
            
            obj.Dofs.(char(name)) = value;
        end
        
        function obj = setCoefficient(obj, name, value)
            % SETCOEFFICIENT Set numeric coefficient data or an evaluation callback.
            %
            %   setCoefficient(obj, name, value) sets a named coefficient
            %   with the specified @a value in the state.
            
            arguments
                obj mdemat.data.variable.State
                name {mustBeTextScalar}
                value {mustBeA(value, {'numeric', 'function_handle'})} = []
            end
            
            obj.Coefs.(char(name)) = value;
        end

        function obj = setHistory(obj, name, value)
            % SETHISTORY Set time series data in state history.
            %
            %   obj = setHistory(obj, name, value) sets a named time
            %   series @a value in the state's history.

            arguments
                obj mdemat.data.variable.State
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
                obj mdemat.data.variable.State
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
        
        function record = snapshot(obj, options)
            %SNAPSHOT Produce a versioned value description without mutable caches.
            %   Callbacks are refused unless omitCallbacks=true; their names then
            %   become mandatory explicit inputs to load/fromRecord.
            arguments
                obj mdemat.data.variable.State
                options.omitCallbacks (1, 1) logical = false
            end
            coefficients = obj.Coefs;
            names = fieldnames(coefficients);
            callbacks = cell(0, 1);
            for index = 1:numel(names)
                if isa(coefficients.(names{index}), 'function_handle')
                    if ~options.omitCallbacks
                        error('mdemat:persistence:UnsupportedCallback', ...
                            'Coefficient %s is callable; use omitCallbacks=true explicitly.', names{index});
                    end
                    callbacks{end+1, 1} = names{index}; %#ok<AGROW>
                    coefficients = rmfield(coefficients, names{index});
                end
            end
            if ~all(structfun(@isnumeric, obj.Dofs)) || ...
                    ~all(structfun(@isnumeric, coefficients)) || ...
                    ~all(structfun(@isnumeric, obj.History))
                error('mdemat:persistence:InvalidState', ...
                    'DOFs, coefficient values and histories must contain numeric arrays.');
            end
            record = struct('format', 'mdemat.state', 'schemaVersion', 1, 'kind', 'base', ...
                'dofs', obj.Dofs, 'coefs', coefficients, 'history', obj.History, ...
                'callbackNames', {callbacks}, 'spaces', struct(), 'parameters', struct());
        end

        function obj = save(obj, filename, options)
            %SAVE Write a data-only snapshot; refuse overwrite by default.
            %   Existing opaque MAT state files are not this versioned format.
            arguments
                obj mdemat.data.variable.State
                filename {mustBeTextScalar}
                options.overwrite (1, 1) logical = false
                options.createDirectory (1, 1) logical = false
                options.omitCallbacks (1, 1) logical = false
            end
            record = obj.snapshot(omitCallbacks=options.omitCallbacks);
            mdemat.tools.io.writeRecord(filename, record, overwrite=options.overwrite, ...
                createDirectory=options.createDirectory);
        end
    end
    
    methods (Static)
        function obj = load(filename, options)
            %LOAD Read this fixed state kind without MATLAB object deserialization.
            arguments
                filename {mustBeTextScalar}
                options.callbacks (1, 1) struct = struct()
            end
            record = mdemat.tools.io.readRecord(filename);
            obj = mdemat.data.variable.State.fromRecord(record, callbacks=options.callbacks);
        end

        function obj = fromRecord(record, options)
            %FROMRECORD Reconstruct a base state from a validated plain record.
            arguments
                record (1, 1) struct
                options.callbacks (1, 1) struct = struct()
            end
            record = mdemat.data.variable.validateStateRecord(record, 'base');
            record = mdemat.data.variable.bindStateCallbacks(record, options.callbacks);
            if ~isempty(fieldnames(record.spaces)) || ~isempty(fieldnames(record.parameters))
                error('mdemat:persistence:InvalidState', 'Base state must not contain spaces/parameters.');
            end
            obj = mdemat.data.variable.State();
            obj.restoreRecordData(record);
        end
    end

    methods (Access = protected)
        function restoreRecordData(obj, record)
            % Internal assignment occurs only after schema/callback validation.
            obj.Dofs = record.dofs;
            obj.Coefs = record.coefs;
            obj.History = record.history;
        end
    end
end
