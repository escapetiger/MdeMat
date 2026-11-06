classdef SpatialState < mdemat.data.variable.State
    % SPATIALSTATE State class with spatial discretization.
    %
    %   SpatialState extends State by adding spatial discretization
    %   information for spatial partial differential equations.
    %
    % See also:
    %   mdemat.data.variable.State, mdemat.common.fem.KineticState

    properties
        XDisc   % Spatial discretization
    end
    
    properties (Dependent)
        nXDims % Number of spatial dimensions
    end

    methods
        function obj = SpatialState(xDisc)
            % SPATIALSTATE Construct an instance of SpatialState.
            %
            %   obj = SpatialState(xDisc) creates a spatial solution state
            %   object with the specified spatial discretization @a xDisc.

            arguments
                xDisc mdemat.common.fem.MeshSpace
            end

            obj@mdemat.data.variable.State();
            obj.XDisc = xDisc;
        end

        function n = get.nXDims(obj)
            % GET.nXDims Get number of spatial dimensions.

            n = obj.XDisc.NDims;
        end

        function dofs = historyDofs(obj, packed)
            % HISTORYDOFS Decode accepted evolution data for observation.
            values = reshape(packed, obj.XDisc.NGlobalDofs, []);
            dofs = struct();
            names = fieldnames(obj.Dofs);
            offset = 0;
            for i = 1:numel(names)
                name = names{i};
                count = size(obj.Dofs.(name), 2);
                dofs.(name) = values(:, offset+(1:count));
                offset = offset+count;
            end
            if offset ~= size(values, 2)
                error('mdemat:state:InvalidHistory', ...
                    'Packed history must contain the declared evolution fields.');
            end
        end

        function newObj = refine(obj, nLevels)
            % REFINE Create refined spatial state.
            %
            %   newObj = refine(obj, nLevels) creates a refined spatial
            %   state with @a nLevels refinement levels.

            arguments
                obj mdemat.common.fem.SpatialState
                nLevels {mustBeNonnegative, mustBeInteger}
            end

            newXDisc = obj.XDisc.refine(nLevels);
            newObj = mdemat.common.fem.SpatialState(newXDisc);
        end

        function obj = report(obj)
            % REPORT Generate spatial state information report.
            %
            %   obj = report(obj) displays the names of degrees of freedom
            %   and coefficients, plus mesh resolution information.

            arguments
                obj mdemat.common.fem.SpatialState
            end
            
            report@mdemat.data.variable.State(obj);

            fprintf('[R] Mesh Resolution: %d elements\n', obj.XDisc.Mesh.NElements);
        end

        function record = snapshot(obj, options)
            %SNAPSHOT Include the current spatial representation without caches.
            arguments
                obj mdemat.common.fem.SpatialState
                options.omitCallbacks (1, 1) logical = false
            end
            record = snapshot@mdemat.data.variable.State(obj, ...
                omitCallbacks=options.omitCallbacks);
            record.kind = 'spatial';
            record.spaces.x = mdemat.common.fem.describeSpace(obj.XDisc);
        end

    end

    methods (Static)
        function obj = load(filename, options)
            %LOAD Read a data-only spatial state with explicit callback binding.
            arguments
                filename {mustBeTextScalar}
                options.callbacks (1, 1) struct = struct()
            end
            record = mdemat.tools.io.readRecord(filename);
            obj = mdemat.common.fem.SpatialState.fromRecord(record, callbacks=options.callbacks);
        end

        function obj = fromRecord(record, options)
            arguments
                record (1, 1) struct
                options.callbacks (1, 1) struct = struct()
            end
            record = mdemat.data.variable.validateStateRecord(record, 'spatial');
            record = mdemat.data.variable.bindStateCallbacks(record, options.callbacks);
            if ~isempty(fieldnames(record.parameters))
                error('mdemat:persistence:InvalidState', 'Spatial state has no additional parameters.');
            end
            xDisc = mdemat.common.fem.SpatialState.restoreSpaces(record, false);
            obj = mdemat.common.fem.SpatialState(xDisc);
            obj.restoreRecordData(record);
        end
    end

    methods (Static, Access = protected)
        function [xDisc, vDisc] = restoreSpaces(record, velocity)
            expected = {'x'};
            if velocity, expected{end+1} = 'v'; end
            if ~isequal(sort(fieldnames(record.spaces)), sort(expected(:)))
                error('mdemat:persistence:InvalidState', 'State space fields do not match its kind.');
            end
            xDisc = mdemat.common.fem.restoreSpace(record.spaces.x);
            vDisc = [];
            if velocity, vDisc = mdemat.common.fem.restoreSpace(record.spaces.v); end
            if isa(xDisc, 'mdemat.common.fem.FiniteElementSpace')
                data = struct2cell(record.dofs);
                valid = cellfun(@(value) isempty(value) || ...
                    size(value, 1) == xDisc.NGlobalDofs, data);
                if ~all(valid)
                    error('mdemat:persistence:InvalidState', ...
                        'Spatial DOF rows must follow the reconstructed local-then-element packing.');
                end
            end
        end
    end
end
