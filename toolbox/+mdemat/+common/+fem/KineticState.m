classdef KineticState < mdemat.common.fem.SpatialState
    % KINETICSTATE State class with spatial and velocity discretizations.
    %
    %   KineticState extends SpatialState by adding velocity discretization
    %   information for kinetic equations such as Boltzmann, Vlasov, or
    %   other velocity-dependent partial differential equations.
    %
    % See also:
    %   mdemat.common.fem.SpatialState, mdemat.data.variable.State

    properties
        VDisc % Velocity discretization
    end

    methods
        function obj = KineticState(xDisc, vDisc)
            % KINETICSTATE Construct an instance of KineticState.
            %
            %   obj = KineticState(xDisc, vDisc) creates a kinetic solution
            %   state object with the specified spatial and velocity
            %   discretizations @a xDisc and @a vDisc.

            arguments
                xDisc mdemat.common.fem.MeshSpace
                vDisc{mustBeA(vDisc, { ...
                    'mdemat.common.fem.SpectralSpace', ...
                    'mdemat.common.fem.SumSpace', ...
                    'mdemat.common.fem.ScaledAffineSpace'})}
            end

            obj@mdemat.common.fem.SpatialState(xDisc);

            obj.VDisc = vDisc;
        end

        function newObj = refine(obj, nLevels)
            % REFINE Create refined kinetic state.
            %
            %   newObj = refine(obj, nLevels) creates a refined kinetic
            %   state by refining only the spatial discretization while
            %   keeping the velocity discretization unchanged.

            arguments
                obj mdemat.common.fem.KineticState
                nLevels{mustBeNonnegative, mustBeInteger}
            end

            newXDisc = obj.XDisc.refine(nLevels);

            newObj = mdemat.common.fem.KineticState(newXDisc, obj.VDisc);
        end

        function obj = report(obj)
            % REPORT Generate kinetic state information report.
            %
            %   obj = report(obj) displays the names of degrees of freedom
            %   and coefficients, plus both spatial and velocity mesh
            %   resolution information.

            arguments
                obj mdemat.common.fem.KineticState
            end

            report@mdemat.common.fem.SpatialState(obj);

            fprintf('[R] Velocity Resolution: %d DoFs\n', obj.VDisc.NDofs);
        end

        function record = snapshot(obj, options)
            %SNAPSHOT Preserve both explicit spatial and velocity descriptions.
            arguments
                obj mdemat.common.fem.KineticState
                options.omitCallbacks (1, 1) logical = false
            end
            record = snapshot@mdemat.common.fem.SpatialState(obj, ...
                omitCallbacks=options.omitCallbacks);
            record.kind = 'kinetic';
            record.spaces.v = mdemat.common.fem.describeSpace(obj.VDisc);
        end

    end

    methods (Static)
        function obj = load(filename, options)
            %LOAD Reconstruct a kinetic state without opaque MATLAB objects.
            arguments
                filename {mustBeTextScalar}
                options.callbacks (1, 1) struct = struct()
            end
            record = mdemat.tools.io.readRecord(filename);
            obj = mdemat.common.fem.KineticState.fromRecord(record, callbacks=options.callbacks);
        end

        function obj = fromRecord(record, options)
            arguments
                record (1, 1) struct
                options.callbacks (1, 1) struct = struct()
            end
            record = mdemat.data.variable.validateStateRecord(record, 'kinetic');
            record = mdemat.data.variable.bindStateCallbacks(record, options.callbacks);
            if ~isempty(fieldnames(record.parameters))
                error('mdemat:persistence:InvalidState', 'Kinetic state has no additional parameters.');
            end
            [xDisc, vDisc] = mdemat.common.fem.SpatialState.restoreSpaces(record, true);
            obj = mdemat.common.fem.KineticState(xDisc, vDisc);
            obj.restoreRecordData(record);
        end
    end
end
