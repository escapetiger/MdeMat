classdef SpatialState < physics.state.State
    % SPATIALSTATE State class with spatial discretization.
    %
    %   SpatialState extends State by adding spatial discretization
    %   information for spatial partial differential equations.
    %
    % See also:
    %   physics.state.State, physics.state.KineticState

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
                xDisc approx.space.MeshSpace
            end

            obj@physics.state.State();
            obj.XDisc = xDisc;
        end

        function n = get.nXDims(obj)
            % GET.nXDims Get number of spatial dimensions.

            n = obj.XDisc.NDims;
        end

        function newObj = refine(obj, nLevels)
            % REFINE Create refined spatial state.
            %
            %   newObj = refine(obj, nLevels) creates a refined spatial
            %   state with @a nLevels refinement levels.

            arguments
                obj physics.state.SpatialState
                nLevels {mustBeNonnegative, mustBeInteger}
            end

            newXDisc = obj.XDisc.refine(nLevels);
            newObj = physics.state.SpatialState(newXDisc);
        end

        function obj = report(obj)
            % REPORT Generate spatial state information report.
            %
            %   obj = report(obj) displays the names of degrees of freedom
            %   and coefficients, plus mesh resolution information.

            arguments
                obj physics.state.SpatialState
            end
            
            report@physics.state.State(obj);

            fprintf('[R] Mesh Resolution: %d elements\n', obj.XDisc.Mesh.NElements);
        end

        function obj = save(obj, filename)
            % SAVE Save spatial state to file.
            %
            %   obj = save(obj, filename) saves the complete spatial state
            %   including discretization data to a MAT file with the
            %   specified @a filename.

            arguments
                obj physics.state.SpatialState
                filename {mustBeTextScalar}
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
            state.className = class(obj);
            
            save(filename, 'state', '-v7.3');
            
            fprintf('[M] Spatial state saved to: %s\n', filename);
        end
    end

    methods (Static)
        function obj = load(filename)
            % LOAD Load spatial state from file.
            %
            %   obj = SpatialState.load(filename) loads a spatial state
            %   object from the specified MAT file and returns a new
            %   SpatialState instance.

            arguments
                filename {mustBeTextScalar}
            end

            core.except.assert(exist(filename, 'file') == 2, ...
                'InvalidInput', 'File %s does not exist', filename);

            state = load(filename).state;
            
            obj = physics.state.SpatialState(state.xDisc);
            
            obj.Dofs = state.dofs;
            obj.Coefs = state.coefs;
            obj.History = state.history;
            
            fprintf('[M] Spatial state loaded from: %s\n', filename);
        end
    end
end