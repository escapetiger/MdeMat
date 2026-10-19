classdef KineticState < physics.SpatialState
    % KINETICSTATE State class with spatial and velocity discretizations.
    %
    %   KineticState extends the SpatialState class by adding velocity
    %   discretization information for kinetic equations such as Boltzmann,
    %   Vlasov, or other velocity-dependent partial differential equations.
    %
    % See Also:
    %   physics.SpatialState, physics.State

    properties
        vDisc %< Velocity discretization
    end

    methods
        function obj = KineticState(xDisc, vDisc)
            % KINETICSTATE Constructor for KineticState.
            %
            %   obj = KineticState(xDisc, vDisc) creates a kinetic solution
            %   state object with the specified spatial and velocity
            %   discretizations and initializes empty degrees of freedom
            %   and coefficients.
            %
            % Inputs:
            %   xDisc - Spatial discretization object
            %   vDisc - Velocity discretization object
            %
            % Outputs:
            %   obj - Constructed KineticState object

            %< Input validation (Rule 235)
            core.except.assert(nargin == 2, 'InvalidInput', ...
                'KineticState requires exactly 2 input arguments');

            %< Call parent constructor
            obj@physics.SpatialState(xDisc);
            
            %< Set velocity discretization
            obj.vDisc = vDisc;
        end

        function newObj = refine(obj, nLevels)
            % REFINE Create refined kinetic state.
            %
            %   newObj = refine(obj, nLevels) creates a refined kinetic
            %   state by refining only the spatial discretization while
            %   keeping the velocity discretization unchanged.
            %
            % Inputs:
            %   obj - The KineticState object
            %   nLevels - Number of spatial refinement levels (non-negative integer)
            %
            % Outputs:
            %   newObj - Refined KineticState object

            %< Input validation (Rule 235)
            core.except.assert(isscalar(nLevels) && nLevels >= 0 && ...
                floor(nLevels) == nLevels, 'InvalidInput', ...
                'nLevels must be a non-negative integer');

            %< Refine only spatial discretization
            newXDisc = obj.xDisc.refine(nLevels);
            
            %< Create new kinetic state with refined spatial discretization
            newObj = physics.KineticState(newXDisc, obj.vDisc);
        end

        function obj = report(obj)
            % REPORT Generate kinetic state information report.
            %
            %   obj = report(obj) displays the names of degrees of freedom
            %   and coefficients, plus both spatial and velocity mesh
            %   resolution information.
            %
            % Inputs:
            %   obj - The KineticState object
            %
            % Outputs:
            %   obj - The KineticState object
            
            %< Call parent report method
            report@physics.SpatialState(obj);

            %< Report velocity discretization information
            if isa(obj.vDisc, 'approx.space.MeshSpace')
                fprintf('[R] Velocity Mesh Resolution: %d elements\n', ...
                    obj.vDisc.mesh.nElements);
            else
                fprintf('[R] Velocity Mesh Resolution: (none)\n');
            end
        end

        function obj = save(obj, filename)
            % SAVE Save kinetic state to file.
            %
            %   obj = save(obj, filename) saves the complete kinetic state
            %   including spatial and velocity discretization data to a
            %   MAT file.
            %
            % Inputs:
            %   obj - The KineticState object
            %   filename - Output filename (char or string)
            %
            % Outputs:
            %   obj - The KineticState object

            %< Input validation (Rule 235)
            core.except.assert(ischar(filename) || isstring(filename), ...
                'InvalidInput', 'filename must be char or string');

            %< Save kinetic state data
            state.dofs = obj.dofs;
            state.coefs = obj.coefs;
            state.xDisc = obj.xDisc;
            state.vDisc = obj.vDisc;
            state.className = class(obj);
            
            save(filename, 'state', '-v7.3');
            
            fprintf('[M] Kinetic state saved to: %s\n', filename);
        end
    end

    methods (Static)
        function obj = load(filename)
            % LOAD Load kinetic state from file.
            %
            %   obj = physics.KineticState.load(filename) loads a kinetic
            %   state object from the specified MAT file and returns a new
            %   KineticState instance.
            %
            % Inputs:
            %   filename - Input filename (char or string)
            %
            % Outputs:
            %   obj - New KineticState object with loaded data

            %< Input validation (Rule 235)
            core.except.assert(ischar(filename) || isstring(filename), ...
                'InvalidInput', 'filename must be char or string');
            core.except.assert(exist(filename, 'file') == 2, ...
                'InvalidInput', 'File %s does not exist', filename);

            %< Load state data
            state = load(filename).state;
            
            %< Validate loaded data has required fields
            requiredFields = {'xDisc', 'vDisc'};
            [tf, mf] = core.validate.hasFields(state, requiredFields);
            core.except.assert(tf, 'InvalidInput', ...
                'Loaded file missing required field: %s', mf);
            
            %< Create new kinetic state object
            obj = physics.KineticState(state.xDisc, state.vDisc);
            
            %< Restore state data
            obj.dofs = state.dofs;
            obj.coefs = state.coefs;
            
            fprintf('[M] Kinetic state loaded from: %s\n', filename);
        end
    end
end