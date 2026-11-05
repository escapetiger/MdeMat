classdef Visualizer < handle
    % VISUALIZER Base class for visualizers.
    %
    %   The Visualizer class provides comprehensive plotting capabilities
    %   for advection equation solutions, including support for exact
    %   solution comparisons and temporal interpolation. It manages plot
    %   timelines, rendering strategies, and coordinate visualization
    %   workflows.
    %
    %   Key features include:
    %   - Multi-dimensional visualization support (1D, 2D, 3D)
    %   - Exact solution overlay capabilities
    %   - Temporal interpolation for smooth animations
    %   - Flexible component-based plotting configuration
    %   - Automatic subplot arrangement and default styling
    %   - Automatic render type selection based on data type
    %
    % Examples:
    %   % Basic 1D visualization setup
    %   visualizer = Visualizer('nDims', 1, 'density', [50], 'final', 1.0);
    %
    %   % With exact solution comparison
    %   visualizer = Visualizer('nDims', 2, 'density', [10, 10], 'final', 2.0, ...
    %                          'components', struct('u', [1]));
    %   visualizer.addExact('u', @(x,y,t) exp(-(x-t).^2-(y-t).^2));
    %   visualizer.plot(space, dofs, tDisc);
    %
    % Notes:
    %   The class uses a handle-based design to manage visualization state
    %   efficiently across multiple time steps and solution updates.
    %   Default styling is automatically applied based on dataset types.
    %
    % See Also:
    %   profilers.visual.Strategy, profilers.visual.GeometryPlotter

    properties
        experimentId % Experiment identifier (string)
        schemeId % Scheme identifier (string)
        nDims % Number of spatial dimensions (scalar)
        timeline % Plot timeline object (approx.mesh.StaticTimeline)
        density % Number of plot points per element (vector)
        components % Plot components structure
        exacts % Exact solution functions (struct)
        plotter % Plotting object (profilers.visual.GeometryPlotter)
        dataset % Visualization dataset (struct)
        isPrepared = false % Flag indicating preparation status (logical)
    end

    properties (Dependent)
        isEnabled % Flag indicating whether visualization is enabled (logical)
        hasExact % Flag indicating whether has exact solution (logical)
        nComponents % Number of components per field (vector)
        nAxesPerComponent % Number of axes per component (scalar)
        nRows % Number of subplot rows (scalar)
        nCols % Number of subplot columns (scalar)
        nAxes % Total number of axes (scalar)
    end

    methods
        function obj = Visualizer(varargin)
            % VISUALIZER Constructor for Visualizer.
            %
            %   obj = Visualizer() creates a default visualizer with empty
            %   configuration requiring manual setup.
            %
            %   obj = Visualizer('Parameter', Value, ...) creates a
            %   visualizer with specified parameters for dimension,
            %   timeline, and plotting configuration.
            %
            % Inputs:
            %   varargin - Input arguments
            %   'experimentId' - Experiment identifier (string, default: '')
            %   'schemeId' - Scheme identifier (string, default: '')
            %   'nDims' - Number of spatial dimensions (positive integer, default: [])
            %   'nTimeNodes' - Number of time nodes for output (positive integer, default: 1)
            %   'final' - Final time for visualization timeline (non-negative scalar, default: 0)
            %   'density' - Plot point density per element (positive vector, default: [])
            %   'components' - Component configuration structure (struct, default: [])
            %
            % Outputs:
            %   obj - Constructed Visualizer object

            p = inputParser;
            addParameter(p, 'experimentId', '', @(x) ischar(x) || isstring(x));
            addParameter(p, 'schemeId', '', @(x) ischar(x) || isstring(x));
            addParameter(p, 'nDims', [], @(x) isnumeric(x) && isscalar(x));
            addParameter(p, 'nTimeNodes', 1, @(x) isnumeric(x) && isscalar(x));
            addParameter(p, 'final', 0, @(x) isnumeric(x) && isscalar(x));
            addParameter(p, 'density', [], @(x) isnumeric(x));
            addParameter(p, 'components', [], @(x) isstruct(x));
            parse(p, varargin{:});

            obj.experimentId = p.Results.experimentId;
            obj.schemeId = p.Results.schemeId;
            obj.nDims = p.Results.nDims;
            obj.timeline = approx.mesh.StaticTimeline(p.Results.nTimeNodes, p.Results.final);
            obj.density = p.Results.density;
            obj.components = p.Results.components;
            obj.exacts = struct();
            obj.dataset = struct();
        end

        function obj = reset(obj)
            % RESET Reset visualizer to initial state.
            %
            %   obj = reset(obj) reinitializes the visualizer by clearing
            %   datasets, resetting the timeline, and reconstructing the
            %   plotting strategy based on current component configuration.
            %   Uses default styling automatically.
            %
            % Inputs:
            %   obj - The Visualizer object
            %
            % Outputs:
            %   obj - The Visualizer object

            strategy = profilers.visual.Strategy(obj.nDims, 'map');
            dofNames = fieldnames(obj.components);
            offset = 0;
            for iField = 1:length(obj.nComponents)
                nFieldComponents = obj.nComponents(iField);
                dofName = dofNames{iField};
                for iComponent = 1:nFieldComponents
                    componentName = sprintf('%s%d', dofName, ...
                        obj.components.(dofName)(iComponent));
                    axisIdx = (offset + iComponent - 1) * obj.nAxesPerComponent + 1;
                    strategy.bind(componentName, axisIdx);
                end
                offset = offset + nFieldComponents;
            end
            obj.plotter = profilers.visual.GeometryPlotter( ...
                core.geometry.Orthotope.unit(obj.nDims), ...
                obj.nRows, obj.nCols);
            obj.plotter.setNodes(obj.density);
            obj.plotter.setStrategy(strategy);
            for iAxis = 1:obj.nAxes
                obj.plotter.addAxisSpec(1, 1);
            end
            obj.plotter.draft();
            obj.timeline.reset();
            obj.isPrepared = false;
        end

        function obj = addDataset(obj, name, type, dataset)
            % ADDDATASET Add dataset for visualization.
            %
            %   obj = addDataset(obj, name, type) registers a new
            %   dataset with the specified name and data type.
            %
            %   obj = addDataset(obj, name, type, dataset) registers a new
            %   dataset with the specified name and data.
            %
            % Inputs:
            %   obj - The Visualizer object
            %   name - Dataset name (string or char array)
            %   type - Dataset type (scalar): 0=fixed, 1=exact, 2=numeric
            %   dataset - Dataset (optional, default: [])
            %
            % Outputs:
            %   obj - The Visualizer object
            %
            % Notes:
            %   Dataset type determines rendering style:
            %   - type < 2: rendered as lines
            %   - type >= 2: rendered as scatter plots

            if nargin < 4 || isempty(dataset)
                obj.dataset.(name) = core.data.GridStruct(type);
                return;
            end
            obj.dataset.(name) = dataset;
            obj.dataset.(name).type = type;
        end

        function obj = addExact(obj, name, exact)
            % ADDEXACT Add exact solution.
            %
            %   obj = addExact(obj, name, exact) registers a new exact
            %   solution with the specified function handle or dataset.
            %
            % Inputs:
            %   obj - The Visualizer object
            %   name - Solution name (string or char array)
            %   exact - Exact function/dataset (function handle or GridStruct)
            %
            % Outputs:
            %   obj - The Visualizer object            
            
            obj.exacts.(name) = exact;
        end

        function obj = plot(obj, space, dofs, tDisc)
            % PLOT Generate visualization for current solution state.
            %
            %   obj = plot(obj, space, dofs, tDisc) creates plots for the
            %   current solution state, including exact solution
            %   comparisons if configured. Handles temporal interpolation
            %   and multi-component visualization automatically using
            %   default styling based on dataset types.
            %
            % Inputs:
            %   obj - The Visualizer object
            %   space - Spatial discretization object
            %   dofs - Degrees of freedom structure containing solution data
            %   tDisc - Temporal discretization object
            %
            % Outputs:
            %   obj - The Visualizer object
            %
            % Notes:
            %   The method automatically handles coordinate preparation on
            %   first call and manages timeline advancement for animations.
            %   Rendering styles are determined by dataset types.

            if ~obj.isPrepared
                x = space.mesh.collocate(obj.plotter.coords);
                dsNames = fieldnames(obj.dataset);
                for i = 1:length(dsNames)
                    dsName = dsNames{i};
                    if obj.dataset.(dsName).type > 0
                        obj.dataset.(dsName).setCoordinates(x);
                    end
                end
                obj.isPrepared = true;
            end

            xRef = obj.plotter.nodes;
            t0 = tDisc.timeline.now;
            t1 = tDisc.timeline.next;
            t = obj.timeline.now;

            dofNames = fieldnames(dofs);

            if t1 > t0
                k = min(tDisc.timeline.nSteps, tDisc.timeline.count);
                V0 = cell(1, k);
                for i = 1:k
                    W0 = reshape(tDisc.U0{i}, space.nGlobalDofs, []);
                    offset = 0;
                    for j = 1:length(dofNames)
                        dofName = dofNames{j};
                        nComponentsPerDof = size(dofs.(dofName), 2);
                        dof = W0(:, offset + (1:nComponentsPerDof));
                        V0{i}.(dofName) = space.evaluate([], xRef, dof);
                        offset = offset + nComponentsPerDof;
                    end
                end
            else
                V0 = [];
            end

            while t1 >= t - 1e-14
                for i = 1:length(dofNames)
                    dofName = dofNames{i};
                    component = obj.components.(dofName);
                    if ~isempty(component)
                        dsNames = fieldnames(obj.dataset);
                        for j = 1:length(dsNames)
                            dsName = dsNames{j};
                            if obj.dataset.(dsName).type == 0
                                continue;
                            end

                            if obj.dataset.(dsName).type == 1
                                f = obj.exacts.(dofName);
                                U = space.functionEvaluate([], xRef, f, t);
                                U = U(:, component);
                            else
                                U = space.evaluate([], xRef, dofs.(dofName)(:, component));
                                if ~isempty(V0)
                                    U0 = arrayfun(@(i) V0{i}.(dofName)(:, component), 1:k, 'Un', 0);
                                    U0 = cellfun(@(x) reshape(x, [1, size(x)]), U0, 'Un', 0);
                                    n = size(U);
                                    U = reshape(U, [1, n]);
                                    U = cat(1, U0{:}, U);
                                    ti = [t1 - cumsum(tDisc.timeline.h(k:-1:1)), t1].';
                                    U = interp1(ti, U, t);
                                    U = reshape(squeeze(U), n);
                                end
                            end
                            obj.dataset.(dsName).data.(dofName) = U;
                        end
                    end
                end

                %< Create title and render
                if isempty(obj.experimentId) && isempty(obj.schemeId)
                    titleText = sprintf('t = %0.2f', t);
                elseif ~isempty(obj.experimentId) && isempty(obj.schemeId)
                    titleText = sprintf('%s at t = %0.2f', obj.experimentId, t);
                elseif ~isempty(obj.schemeId) && isempty(obj.experimentId)
                    titleText = sprintf('%s at t = %0.2f', obj.schemeId, t);
                else
                    titleText = sprintf('%s with %s at t = %0.2f', obj.experimentId, obj.schemeId, t);
                end

                %< Render with default styling
                obj.plotter.render(titleText, obj.dataset);

                %< Advance output timeline
                obj.timeline.advance();
                t = obj.timeline.now;
            end
        end

        function TF = get.isEnabled(obj)
            % GET.ISENABLED Determine if visualization is enabled.

            TF = any(obj.density > 0) && sum(obj.nComponents) > 0;
        end

        function TF = get.hasExact(obj)
            % GET.HASEXACT Determine if exact solutions are available.

            TF = ~isempty(fieldnames(obj.exacts));
        end

        function n = get.nComponents(obj)
            % GET.NCOMPONENTS Get number of components per field.

            if isempty(obj.components)
                n = [];
                return;
            end
            n = cellfun(@(f) length(obj.components.(f)), fieldnames(obj.components));
        end

        function n = get.nAxesPerComponent(obj)
            % GET.NAXESPERCOMPONENT Get number of axes per component.

            switch obj.nDims
                case 1
                    n = 1;
                case {2, 3}
                    n = cellfun(@(f) ~isempty(obj.dataset.(f)), fieldnames(obj.dataset));
                    n = sum(n) + 1;
            end
        end

        function n = get.nRows(obj)
            % GET.NROWS Get number of subplot rows.

            if isempty(obj.nComponents)
                n = 0;
                return;
            end
            
            switch obj.nDims
                case 1
                    n = sum(arrayfun(@(x) min(1, x), obj.nComponents));
                case {2, 3}
                    n = sum(obj.nComponents);
            end
        end

        function n = get.nCols(obj)
            % GET.NCOLS Get number of subplot columns.

            if isempty(obj.nComponents)
                n = 0;
                return;
            end
            
            switch obj.nDims
                case 1
                    n = max(obj.nComponents, [], 1);
                case {2, 3}
                    if isempty(obj.dataset)
                        n = 1;
                    else
                        n = cellfun(@(f) ~isempty(obj.dataset.(f)), fieldnames(obj.dataset));
                        n = sum(n) + 1;
                    end
            end
        end

        function n = get.nAxes(obj)
            % GET.NAXES Get total number of axes.

            if isempty(obj.nComponents)
                n = 0;
                return;
            end
            n = sum(obj.nComponents) * obj.nAxesPerComponent;
        end
    end
end