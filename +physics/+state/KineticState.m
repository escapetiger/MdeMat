classdef KineticState < physics.state.SpatialState
    % KINETICSTATE State class with spatial and velocity discretizations.
    %
    %   KineticState extends SpatialState by adding velocity discretization
    %   information for kinetic equations such as Boltzmann, Vlasov, or
    %   other velocity-dependent partial differential equations.
    %
    % See also:
    %   physics.state.SpatialState, physics.state.State

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
                xDisc approx.space.MeshSpace
                vDisc{mustBeA(vDisc, { ...
                    'approx.space.SpectralSpace', ...
                    'approx.space.SumSpace', ...
                    'approx.space.ScaledAffineSpace'})}
            end

            obj@physics.state.SpatialState(xDisc);

            obj.VDisc = vDisc;
        end

        function newObj = refine(obj, nLevels)
            % REFINE Create refined kinetic state.
            %
            %   newObj = refine(obj, nLevels) creates a refined kinetic
            %   state by refining only the spatial discretization while
            %   keeping the velocity discretization unchanged.

            arguments
                obj physics.state.KineticState
                nLevels{mustBeNonnegative, mustBeInteger}
            end

            newXDisc = obj.XDisc.refine(nLevels);

            newObj = physics.state.KineticState(newXDisc, obj.VDisc);
        end

        function obj = report(obj)
            % REPORT Generate kinetic state information report.
            %
            %   obj = report(obj) displays the names of degrees of freedom
            %   and coefficients, plus both spatial and velocity mesh
            %   resolution information.

            arguments
                obj physics.state.KineticState
            end

            report@physics.state.SpatialState(obj);

            fprintf('[R] Velocity Resolution: %d DoFs\n', obj.VDisc.NDofs);
        end

        function obj = save(obj, filename)
            % SAVE Save kinetic state to file.
            %
            %   obj = save(obj, filename) saves the complete kinetic state
            %   including spatial and velocity discretization data to a MAT
            %   file with the specified @a filename.

            arguments
                obj physics.state.KineticState
                filename{mustBeTextScalar}
            end

            [dirPath, ~, ~] = fileparts(filename);
            if ~isempty(dirPath) && ~exist(dirPath, 'dir')
                mkdir(dirPath);
                fprintf('[M] Created directory: %s\n', dirPath);
            end

            state.dofs = obj.Dofs;
            state.coefs = obj.Coefs;
            state.history = obj.History;
            state.xDisc = obj.XDisc;
            state.vDisc = obj.VDisc;
            state.className = class(obj);

            save(filename, 'state', '-v7.3');

            fprintf('[M] Kinetic state saved to: %s\n', filename);
        end
    end

    methods (Static)
        function obj = load(filename)
            % LOAD Load kinetic state from file.
            %
            %   obj = KineticState.load(filename) loads a kinetic state
            %   object from the specified MAT file and returns a new
            %   KineticState instance.

            arguments
                filename{mustBeTextScalar}
            end

            core.except.assert(exist(filename, 'file') == 2, ...
                'InvalidInput', 'File %s does not exist', filename);

            state = load(filename).state;

            requiredFields = {'xDisc', 'vDisc'};
            [tf, mf] = core.except.hasFields(state, requiredFields);
            core.except.assert(tf, 'InvalidInput', ...
                'Loaded file missing required field: %s', mf);

            obj = physics.state.KineticState(state.xDisc, state.vDisc);

            obj.Dofs = state.dofs;
            obj.Coefs = state.coefs;
            obj.History = state.history;

            fprintf('[M] Kinetic state loaded from: %s\n', filename);
        end
    end
end