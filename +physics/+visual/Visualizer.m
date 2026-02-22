classdef Visualizer < handle
    % VISUALIZER Base class for visualizers.
    %
    %   Visualizer provides comprehensive plotting capabilities for physics
    %   simulations, including support for exact solution comparisons and
    %   temporal interpolation. It manages plot timelines, rendering
    %   strategies, and coordinate visualization workflows using
    %   FigureManager for automatic figure positioning and Database for
    %   centralized data and geometry management.
    %
    %   The class uses a handle-based design to manage visualization state
    %   efficiently across multiple time steps and solution updates.
    %   Default styling is automatically applied based on dataset types.
    %
    % See also:
    %   physics.visual.Strategy, physics.visual.Database,
    %   physics.visual.FigureManager, physics.visual.TiledPlotter

    properties (Constant)
        IsInitialized = 0b100 % Status mask: initialized
        IsRendering = 0b010 % Status mask: plotting
        IsFinalized = 0b001 % Status mask: finalized
    end

    properties (SetAccess = protected, GetAccess = public)
        FigureManager % Figure manager for automatic positioning (FigureManager)
        Database % Centralized data management (Database)
        Plotters % Plotting objects mapped to figure handles (struct)
        Timeline % Plot timeline object (approx.mesh.StaticTimeline)
        Density % Number of plot points per element (vector)
        Components % Plot components structure
        Exacts % Exact solution functions (struct)
        TitlePrefix % Plot title head (string)
        Status % Binary mask indicating the status of visualizer
        Cache % Cache for reference nodes and coordinates (struct)
    end

    properties (Dependent)
        IsEnabled % Flag indicating whether visualization is enabled (logical)
        HasExact % Flag indicating whether has exact solution (logical)
        NComponents % Number of components per field (vector)
    end

    methods
        function obj = Visualizer()
            % VISUALIZER Construct an instance of Visualizer.
            %
            %   obj = Visualizer() creates a Visualizer with default settings.
            obj.Timeline = [];
            obj.Density = 0;
            obj.Components = struct();
            obj.Exacts = struct();
            obj.Plotters = struct();
            obj.TitlePrefix = '';
            obj.Status = 0b000;
            obj.Cache = struct();
            obj.FigureManager = physics.visual.FigureManager();
            obj.Database = physics.visual.Database();
        end

        function obj = setDensity(obj, density)
            % SETDENSITY Set the number of points per element to plot.
            %
            %   obj = setDensity(obj, density) sets the number of points
            %   per element for plotting visualization.

            arguments
                obj physics.visual.Visualizer
                density {mustBeNonnegative, mustBeVector}
            end

            obj.Density = density;
        end

        function obj = setRefXNodes(obj, geometry)
            % SETREFNODES Configure spatial evaluation nodes for rendering.
            %
            %   obj = setRefXNodes(obj, geometry) sets up a uniform grid of
            %   evaluation nodes within the reference geometry bounds using
            %   the provided @a geometry.

            arguments
                obj physics.visual.Visualizer
                geometry core.geometry.Geometry
            end

            if isa(geometry, 'core.geometry.Orthotope')
                nd = geometry.NDims;

                n = floor(obj.Density^(1/nd));
                
                if isscalar(n)
                    n = repmat(n, 1, nd);
                end
                
                a = geometry.Lower;
                b = geometry.Upper;
                h = (b - a) ./ (n + 1);
                x = arrayfun(@(i) (a(i)+h(i)):h(i):(b(i)-h(i)), 1:nd, 'Un', false);
                obj.Cache.RefXCoords = x;

                [X{1:nd}] = ndgrid(x{:});
                obj.Cache.RefXNodes = reshape(cat(nd+1, X{:}), [], nd).';
            else
                core.except.assert(0, 'InvalidGeometry', ...
                    'Unsupported geometry type: %s', class(geometry));
            end
        end

        function obj = setPhyXNodes(obj, mesh)
            % SETPHYNODES Set up physical spatial nodes for the visualizer.
            %
            %   obj = setPhyXNodes(obj, mesh) prepares the visualizer with the
            %   provided mesh by generating physical coordinates and
            %   initializing datasets.

            arguments
                obj physics.visual.Visualizer
                mesh approx.mesh.Mesh
            end

            datasetNames = obj.Database.DatasetNames;
            for j = 1:obj.Database.NDatasets
                datasetName = datasetNames{j};
                dataset = obj.Database.getDataset(datasetName);

                if strcmp(dataset.Type, 'fixed')
                    continue;
                end

                x = mesh.collocate(obj.Cache.RefXNodes);
                dataset.setData('x', x);
            end
        end

        function obj = setPhyVNodes(obj, nodes)
            % SETPHYVNODES Set up physical velocity nodes for the
            % visualizer.
            %
            %   obj = setPhyVNodes(obj, nodes) prepares the visualizer with
            %   the provided velocity nodes.

            arguments
                obj physics.visual.Visualizer
                nodes {mustBeNumeric}
            end

            datasetNames = obj.Database.DatasetNames;
            for j = 1:obj.Database.NDatasets
                datasetName = datasetNames{j};
                dataset = obj.Database.getDataset(datasetName);

                if strcmp(dataset.Type, 'fixed')
                    continue;
                end

                dataset.setData('v', nodes);
            end
        end

        function obj = setTimeline(obj, final, nTimeNodes)
            % SETTIMELINE Set the plotting timeline.
            %
            %   obj = setTimeline(obj, final, nTimeNodes) sets the plotting
            %   timeline with the specified final time and number of time
            %   nodes.

            arguments
                obj physics.visual.Visualizer
                final {mustBeNonnegative}
                nTimeNodes {mustBePositive, mustBeInteger}
            end

            obj.Timeline = approx.mesh.StaticTimeline(nTimeNodes, final);
        end

        function obj = setComponents(obj, components)
            % SETCOMPONENTS Set the solution components to plot.
            %
            %   obj = setComponents(obj, components) sets the indices of
            %   solution components to plot.

            arguments
                obj physics.visual.Visualizer
                components struct
            end

            obj.Components = components;
        end

        function obj = setTitlePrefix(obj, titlePrefix)
            % SETTITLEHEAD Set the prefix of title.
            %
            %   obj = setTitlePrefix(obj, titlePrefix) sets the formatted
            %   title prefix for plotting.

            arguments
                obj physics.visual.Visualizer
                titlePrefix string
            end

            obj.TitlePrefix = strrep(titlePrefix, '_', '-');
        end

        function obj = addDataset(obj, name, type, data)
            % ADDDATASET Add dataset for visualization.
            %
            %   obj = addDataset(obj, name, type) registers a new dataset
            %   with the specified @a name and data @a type in the
            %   centralized Database.
            %
            %   obj = addDataset(obj, name, type, data) registers a new
            %   dataset with initial @a data.
            %
            %   Dataset type determines rendering style in 1D plot:
            %   - 'fixed'/'exact': rendered as lines
            %   - 'numeric': rendered as scatter plots

            arguments
                obj physics.visual.Visualizer
                name {mustBeTextScalar}
                type {mustBeTextScalar}
                data struct = struct()
            end

            dataset = physics.visual.Dataset(Type=type, Data=data);
            obj.Database.setDataset(name, dataset);
        end

        function obj = addExact(obj, name, func)
            % ADDEXACT Add exact solution.
            %
            %   obj = addExact(obj, name, func) registers a new exact
            %   solution with the specified function handle.

            arguments
                obj physics.visual.Visualizer
                name {mustBeTextScalar}
                func {mustBeA(func, {'numeric', 'function_handle'})}
            end

            obj.Exacts.(char(name)) = func;
        end

        function obj = addPlotter(obj, name, strategy)
            % ADDPLOTTER Add tiled plotter.
            %
            %   obj = addPlotter(obj, name, strategy) registers a plotter
            %   with the specified rendering strategy.

            arguments
                obj physics.visual.Visualizer
                name {mustBeTextScalar}
                strategy {mustBeTextScalar}
            end

            switch lower(strategy)
                case '1d'
                    nRows = sum(arrayfun(@(x) min(1, x), obj.NComponents));
                    nCols = max(obj.NComponents, [], 1);
                    nAxes = sum(obj.NComponents);
                    nAxesPerComponent = 1;
                    strategyObj = physics.visual.Strategy1d();
                case '2d'
                    nRows = sum(obj.NComponents);
                    nCols = length(fieldnames(obj.Database.Datasets));
                    nAxesPerComponent = nCols;
                    nAxes = nRows * nCols;
                    strategyObj = physics.visual.Strategy2d();
                case '3d'
                    nRows = sum(obj.NComponents);
                    nCols = length(fieldnames(obj.Database.Datasets));
                    nAxesPerComponent = nCols;
                    nAxes = nRows * nCols;
                    strategyObj = physics.visual.Strategy3d();
                case 'slice1d'
                    nRows = sum(arrayfun(@(x) min(1, x), obj.NComponents));
                    nCols = max(obj.NComponents, [], 1);
                    nAxes = sum(obj.NComponents);
                    nAxesPerComponent = 1;
                    strategyObj = physics.visual.StrategySlice1d();
                otherwise
                    core.except.assert(false, 'InvalidStrategy', ...
                        'Unknown strategy: %s', strategy);
            end

            plotter = physics.visual.TiledPlotter(nRows, nCols);

            dofNames = fieldnames(obj.Components);
            offset = 0;
            for iField = 1:length(obj.NComponents)
                nFieldComponents = obj.NComponents(iField);
                dofName = dofNames{iField};
                for iComponent = 1:nFieldComponents
                    componentName = sprintf('%s%d', dofName, ...
                        obj.Components.(dofName)(iComponent));
                    axisIdx = (offset + iComponent - 1) * nAxesPerComponent + 1;
                    plotter.bind(componentName, axisIdx);
                end
                offset = offset + nFieldComponents;
            end

            plotter.setStrategy(strategyObj);

            for iAxis = 1:nAxes
                plotter.addAxisSpec(1, 1);
            end

            obj.Plotters.(name) = plotter;
        end

        function obj = initialize(obj)
            % INITIALIZE Set up plotters and figures based on dimensions
            % and components.
            %
            %   obj = initialize(obj) creates the appropriate number of
            %   figures via FigureManager and sets up plotters for the
            %   configured dimensions and components. Also configures the
            %   Database with reference geometry for coordinate generation.

            arguments
                obj physics.visual.Visualizer
            end

            if ~isempty(obj.Timeline)
                obj.Timeline.reset();
            end

            if bitand(obj.Status, obj.IsInitialized), return; end

            if obj.IsEnabled
                plotterNames = fieldnames(obj.Plotters);
                nPlotters = length(plotterNames);
                if isempty(obj.FigureManager.FigureHandles)
                    figureHandles = obj.FigureManager.createFigures(nPlotters);
                else
                    figureHandles = obj.FigureManager.FigureHandles;
                end
    
                for iPlotter = 1:nPlotters
                    plotterName = plotterNames{iPlotter};
                    plotter = obj.Plotters.(plotterName);
                    plotter.setFigure(figureHandles{iPlotter});
                    plotter.draft();
                end
            end

            obj.Status = bitor(obj.Status, obj.IsInitialized);
        end

        function obj = update(obj, space, dofs, options)
            % UPDATE Update datasets with current solution.

            arguments
                obj physics.visual.Visualizer
                space approx.space.FiniteElementSpace
                dofs struct
                options.timeline approx.mesh.DynamicTimeline
                options.V0 = []
                options.n0 = 0
            end

            xRef = obj.Cache.RefXNodes;
            V0 = options.V0;
            n0 = options.n0;

            if ~isempty(obj.Timeline)
                t = obj.Timeline.Now;
            else
                t = [];
            end

            dofNames = fieldnames(dofs);

            for i = 1:length(dofNames)
                dofName = dofNames{i};
                if ~isfield(obj.Components, dofName)
                    continue;
                end

                component = obj.Components.(dofName);
                if isempty(component)
                    continue;
                end

                datasetNames = obj.Database.DatasetNames;
                for j = 1:length(datasetNames)
                    datasetName = datasetNames{j};
                    dataset = obj.Database.getDataset(datasetName);

                    if strcmp(dataset.Type, 'fixed')
                        continue;
                    end

                    if strcmp(dataset.Type, 'exact')
                        if isfield(obj.Exacts, dofName)
                            f = obj.Exacts.(dofName);
                            if ~isempty(t)
                                args = {'t', t};
                            else
                                args = {};
                            end
                            U = space.feval(f, xRef, args=args);
                            U = U(:, component);
                        else
                            continue;
                        end
                    else
                        U = dofs.(dofName)(:, component);
                        U = space.eval(xRef, U);
                        if ~isempty(V0)
                            ht = options.timeline.StepSizeQueue;
                            t1 = options.timeline.Next;
                            U0 = cell(1, n0);
                            for k = 1:n0
                                U0{k} = V0{k}.(dofName)(:, component);
                                U0{k} = space.eval(xRef, U0{k});
                                U0{k} = reshape(U0{k}, [1, size(U0{k})]);
                            end
                            n = size(U);
                            U = reshape(U, [1, n]);
                            U = cat(1, U0{:}, U);
                            ti = [t1 - cumsum(ht(n0:-1:1)), t1].';
                            U = interp1(ti, U, t, 'linear', 'extrap');
                            U = reshape(squeeze(U), n);
                        end
                    end

                    dataset.setData(dofName, U);
                end
            end
        end

        function obj = render(obj)
            % RENDER Render all plotters with shared title.

            obj.Status = bitor(obj.Status, obj.IsRendering);

            if ~isempty(obj.Timeline)
                t = obj.Timeline.Now;
                titleText = sprintf('%s at t = %0.2f', obj.TitlePrefix, t);
            else
                titleText = obj.TitlePrefix;
            end
    
            plotterNames = fieldnames(obj.Plotters);
            for i = 1:length(plotterNames)
                plotterName = plotterNames{i};
                plotter = obj.Plotters.(plotterName);
                plotter.render(obj.Database, titleText);
            end

            obj.Status = bitand(obj.Status, bitcmp(obj.IsRendering));
        end

        function obj = finalize(obj)
            % FINALIZE Clean up all plotters and figures.
            %
            %   obj = finalize(obj) closes all plotters, clears visualization
            %   resources, and closes managed figures.

            arguments
                obj physics.visual.Visualizer
            end

            plotterNames = fieldnames(obj.Plotters);
            for i = 1:length(plotterNames)
                plotterName = plotterNames{i};
                obj.Plotters.(plotterName).close();
            end
            obj.Plotters = struct();

            obj.FigureManager.closeAll();

            obj.Database.Datasets = struct();

            obj.Status = bitor(obj.Status, obj.IsFinalized);
        end

        function TF = get.IsEnabled(obj)
            % GET.ISENABLED Determine if visualization is enabled.

            arguments
                obj physics.visual.Visualizer
            end

            TF = false;
            
            if isempty(obj.Density) || ~isnumeric(obj.Density) || any(obj.Density <= 0)
                return;
            end
            
            if isempty(obj.Components) || ~isstruct(obj.Components)
                return;
            else
                names = fieldnames(obj.Components);
                if isempty(names)
                    return;
                end
                flag = true;
                for i = 1:length(names)
                    name = names{i};
                    if ~isempty(obj.Components.(name))
                        flag = false;
                        break;
                    end
                end
                if flag, return; end
            end

            if bitand(obj.Status, obj.IsInitialized)
                if isempty(obj.Plotters) || ~isstruct(obj.Plotters)
                    return;
                end
                    
                plotterNames = fieldnames(obj.Plotters);
                if isempty(plotterNames)
                    return;
                end
            end
            
            TF = true;
        end

        function TF = get.HasExact(obj)
            % GET.HASEXACT Determine if exact solutions are available.

            TF = ~isempty(fieldnames(obj.Exacts));
        end

        function n = get.NComponents(obj)
            % GET.NCOMPONENTS Get number of components per field.

             if isempty(obj.Components) || ~isstruct(obj.Components)
                n = [];
                return;
            end
            fieldNames = fieldnames(obj.Components);
            n = cellfun(@(f) length(obj.Components.(f)), fieldNames);
        end
    end
end