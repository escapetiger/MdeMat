classdef Strategy < handle
    % STRATEGY Abstract base class for data rendering strategies.
    %
    %   Strategy defines the interface and common functionality for all
    %   rendering strategies used in data visualization. This class
    %   provides mode-based rendering control, variable-to-axis mapping
    %   capabilities, and utility methods for coordinated multi-panel
    %   visualizations.
    %
    %   Concrete strategy implementations handle the specific details of
    %   rendering data in 1D, 2D, or 3D contexts, while this base class
    %   manages the overall rendering workflow and variable organization.
    %
    % See also:
    %   profilers.visual.Strategy1d, profilers.visual.Strategy2d,
    %   profilers.visual.Strategy3d, profilers.visual.TiledPlotter

    properties
        nDims   % Number of spatial dimensions handled by this strategy
        varMap  % Structure mapping variable names to axis indices
        mode    % Rendering mode ('sequential' | 'map')
    end

    methods
        function obj = Strategy(nDims, mode)
            % STRATEGY Constructor for Strategy base class.
            %
            %   obj = Strategy(nDims, mode) creates a rendering strategy
            %   for the specified number of dimensions with the given
            %   rendering mode.
            %
            % Inputs:
            %   nDims - Number of spatial dimensions (1, 2, or 3)
            %   mode - Rendering mode
            %
            % Notes:
            %   mode = 'sequential': Variables rendered in order of appearance
            %   mode = 'map': Variables rendered according to bind() mapping
            %
            % Outputs:
            %   obj - Constructed Strategy object

            if nargin < 2, mode = 'sequential'; end

            core.except.assert(ismember(mode, {'sequential', 'map'}), ...
                'InvalidInput', 'Mode must be ''sequential'' or ''map''');

            obj.nDims = nDims;
            obj.varMap = struct();
            obj.mode = mode;
        end

        function obj = bind(obj, varName, axisIdx)
            % BIND Map a variable to a specific axis index.
            %
            %   obj = bind(obj, varName, axisIdx) creates a mapping from
            %   the specified variable name to an axis index for use in
            %   'map' rendering mode. This allows precise control over
            %   which variables appear on which axes.
            %
            % Inputs:
            %   obj - The Strategy object
            %   varName - Name of the variable to map (string)
            %   axisIdx - Target axis index (positive integer)
            % 
            % Outputs:
            %   obj - The Strategy object

            obj.varMap.(varName) = axisIdx;
        end

        function obj = render(obj, axes, dataset, styleset)
            % RENDER Render data with variable-to-axis coordination.
            %
            %   obj = render(obj, axes, dataset, styleset) renders the
            %   provided dataset onto the axes array using either
            %   sequential or mapped variable placement based on the
            %   strategy mode.
            %
            % Inputs:
            %   obj - The Strategy object
            %   axes - Array of axes handles for rendering
            %   dataset - Structure containing data fields to render
            %   styleset - Structure containing rendering style specifications
            % 
            % Outputs:
            %   obj - The Strategy object

            if strcmp(obj.mode, 'sequential')
                varNames = fieldnames(dataset);
                obj.renderSequential(axes, dataset, styleset, varNames);
            else
                varNames = fieldnames(obj.varMap);
                obj.renderMapped(axes, dataset, styleset, varNames);
            end
        end
    end

    methods (Access = private)
        function obj = renderSequential(obj, axes, dataset, styleset, varNames)
            % RENDERSEQUENTIAL Render variables in sequential order.

            for j = 1:length(varNames)
                varName = varNames{j};

                nAxes = length(fieldnames(dataset)) + 1;
                offset = 1 + (j - 1) * nAxes;

                if offset + nAxes - 1 > length(axes)
                    core.except.verify(false, 'NotEnoughAxes', ...
                        'Not enough axes for all variables.');
                    break;
                end

                obj.renderVariable(axes, dataset, styleset, varName, offset);
            end
        end

        function renderMapped(obj, axes, dataset, styleset, varNames)
            % RENDERMAPPED Render variables using explicit axis mapping.

            for j = 1:length(varNames)
                varName = varNames{j};
                offset = obj.varMap.(varName);

                if offset > length(axes)
                    core.except.verify(false, 'InvalidAxisIdx', ...
                        'Axis index %d exceeds available axes.', offset);
                    continue;
                end

                obj.renderVariable(axes, dataset, styleset, varName, offset);
            end
        end
    
        function renderVariable(obj, axes, dataset, styleset, varName, offset)
            % RENDERVARIABLE Render a single variable onto axes.

            switch obj.nDims
                case 1
                    renderVariable1d(axes, dataset, styleset, varName, offset);
                case 2
                    renderVariable2d(axes, dataset, styleset, varName, offset);
                case 3
                    renderVariable3d(axes, dataset, styleset, varName, offset);
            end
        end

    end
end