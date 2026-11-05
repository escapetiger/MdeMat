classdef SpatialState < physics.State
    % SPATIALSTATE State class with spatial discretization.
    %
    %   SpatialState extends the base State class by adding spatial
    %   discretization information.
    %
    % See Also:
    %   physics.State

    properties
        xDisc %< Spatial discretization
    end

    methods
        function obj = SpatialState(xDisc)
            % SPATIALSTATE Constructor for SpatialState.
            %
            %   obj = SpatialState(xDisc) creates a spatial solution state
            %   object with the specified spatial discretization and
            %   initializes empty degrees of freedom and coefficients.
            %
            % Inputs:
            %   xDisc - Spatial discretization object
            %
            % Outputs:
            %   obj - Constructed SpatialState object

            obj@physics.State();
            obj.xDisc = xDisc;
        end

        function newObj = refine(obj, nLevels)
            % REFINE Create refined spatial state.
            %
            %   newObj = refine(obj, nLevels) creates a refined spatial
            %   state.
            %
            % Inputs:
            %   obj - The SpatialState object
            %   nLevels - Number of refinement levels (non-negative integer)
            %
            % Outputs:
            %   newObj - Refined SpatialState object

            newXDisc = obj.xDisc.refine(nLevels);
            newObj = physics.SpatialState(newXDisc);
        end

        function obj = report(obj)
            % REPORT Generate spatial state information report.
            %
            %   obj = report(obj) displays the names of degrees of freedom
            %   and coefficients, plus mesh resolution information.
            %
            % Inputs:
            %   obj - The SpatialState object
            %
            % Outputs:
            %   obj - The SpatialState object
            
            report@physics.State(obj);

            if isa(obj.xDisc, 'approx.space.MeshSpace')
                fprintf('[R] Mesh Resolution: %d elements\n', ...
                    obj.xDisc.mesh.nElements);
            else
                fprintf('[R] Mesh Resolution: (none)\n');
            end
        end

        function obj = save(obj, filename)
            % SAVE Save spatial state to file.
            %
            %   obj = save(obj, filename) saves the complete spatial state
            %   including discretization data to a MAT file.
            %
            % Inputs:
            %   obj - The SpatialState object
            %   filename - Output filename (char or string)
            %
            % Outputs:
            %   obj - The SpatialState object

            %< Input validation (rule 235)
            core.except.assert(ischar(filename) || isstring(filename), ...
                'InvalidInput', 'filename must be char or string');

            %< Save spatial state data
            state.dofs = obj.dofs;
            state.coefs = obj.coefs;
            state.xDisc = obj.xDisc;
            state.className = class(obj);
            
            save(filename, 'state', '-v7.3');
            
            fprintf('[M] Spatial state saved to: %s\n', filename);
        end
    end

    methods (Static)
        function obj = load(filename)
            % LOAD Load spatial state from file.
            %
            %   obj = physics.SpatialState.load(filename) loads a spatial
            %   state object from the specified MAT file and returns a new
            %   SpatialState instance.
            %
            % Inputs:
            %   filename - Input filename (char or string)
            %
            % Outputs:
            %   obj - New SpatialState object with loaded data

            %< Input validation (rule 235)
            core.except.assert(ischar(filename) || isstring(filename), ...
                'InvalidInput', 'filename must be char or string');
            core.except.assert(exist(filename, 'file') == 2, ...
                'InvalidInput', 'File %s does not exist', filename);

            %< Load state data
            state = load(filename).state;
            
            %< Create new spatial state object
            obj = physics.SpatialState(state.xDisc);
            
            %< Restore state data
            obj.dofs = state.dofs;
            obj.coefs = state.coefs;
            
            fprintf('[M] Spatial state loaded from: %s\n', filename);
        end
    end
end