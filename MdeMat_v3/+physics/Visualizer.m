classdef Visualizer < handle
    % VISUALIZER Base class for visualizers.
    %
    %   The Visualizer class provides comprehensive plotting capabilities,
    %   including support for exact solution comparisons and temporal
    %   interpolation. It manages plot timelines, rendering strategies, and
    %   coordinate visualization workflows using FigureManager for
    %   automatic figure positioning and DataSource for centralized data
    %   and geometry management.
    %
    % Notes:
    %   The class uses a handle-based design to manage visualization state
    %   efficiently across multiple time steps and solution updates.
    %   Default styling is automatically applied based on dataset types.
    %
    % See Also:
    %   profilers.visual.Strategy, profilers.visual.Database,
    %   profilers.visual.FigureManager, profilers.visual.TiledPlotter

    properties (Constant)
        INITIALIZED = 0b100 % Status mask: initialized
        PLOTTING = 0b010 % Status mask: plotting
        FINALIZED = 0b001 % Status mask: finalized
    end

    properties (SetAccess = protected, GetAccess = public)
        figureManager % Figure manager for automatic positioning (FigureManager)
        database % Centralized data management (Database)
        plotters % Plotting objects mapped to figure handles (struct)
        timeline % Plot timeline object (approx.mesh.StaticTimeline)
        nDims % Number of spatial dimensions (scalar)
        density % Number of plot points per element (vector)
        components % Plot components structure
        exacts % Exact solution functions (struct)
        titleHead % Plot title head (string)
        status % Binary mask indicating the status of visualizer
    end

    properties (Dependent)
        isEnabled % Flag indicating whether visualization is enabled (logical)
        hasExact % Flag indicating whether has exact solution (logical)
        nComponents % Number of components per field (vector)
    end

    methods
        function obj = Visualizer(nDims)
            % VISUALIZER Constructor for Visualizer.

            obj.nDims = nDims;
            obj.timeline = [];
            obj.density = [];
            obj.components = struct();
            obj.exacts = struct();
            obj.plotters = struct();
            obj.titleHead = '';
            obj.status = 0b000;
            obj.figureManager = profilers.visual.FigureManager();
            obj.database = profilers.visual.Database();
        end

        function obj = setDensity(obj, density)
            % SETDENSITY Set the number of points per element to plot.
            %
            %   obj = setDensity(obj, density) sets the number of points
            %   per element for plotting.
            %
            % Inputs:
            %   obj - The Visualizer object
            %   density - Number of points per elements (vector)
            %
            % Outputs:
            %   obj - The Visualizer object

            obj.density = density;
        end

        function obj = setTimeline(obj, final, nTimeNodes)
            % SETTIMELINE Set the plotting timeline.
            %
            %   obj = setTimeline(obj, final) sets the plotting timeline
            %   with the specified final time and number of time nodes.
            %
            % Inputs:
            %   obj - The Visualizer object
            %   final - Final time (scalar)
            %   nTimeNodes - Number of time nodes to plot (scalar)
            %
            % Outputs:
            %   obj - The Visualizer object

            obj.timeline = approx.mesh.StaticTimeline(nTimeNodes, final);
        end

        function obj = setComponents(obj, components)
            % SETCOMPONENTS Set the solution components to plot.
            %
            %   obj = setComponents(obj, components) sets the indices of
            %   solution components to plot.
            %
            % Inputs:
            %   obj - The Visualizer object
            %   components - Indices of solution components (struct)
            %
            % Outputs:
            %   obj - The Visualizer object

            obj.components = components;
        end

        function obj = setTitleHead(obj, titleHead)
            % SETTITLEHEAD Add the head of title.
            %
            %   obj = setTitleHead(obj) sets the formatted title head
            %   string for plotting.
            %
            % Inputs:
            %   obj - The Visualizer object
            %   titleHead - Title head (string)
            %
            % Outputs:
            %   obj - The Visualizer object

            obj.titleHead = strrep(titleHead, '_', '-');
        end

        function obj = addDataset(obj, name, type, data)
            % ADDDATASET Add dataset for visualization.
            %
            %   obj = addDataset(obj, name, type) registers a new
            %   dataset with the specified name and data type in the
            %   centralized Database.
            %
            % Inputs:
            %   obj - The Visualizer object
            %   name - Dataset name (string or char array)
            %   type - Dataset type (scalar): 0=fixed, 1=exact, 2=numeric
            %   data - Data structure (struct, optional)
            %
            % Outputs:
            %   obj - The Visualizer object
            %
            % Notes:
            %   Dataset type determines rendering style:
            %   - type < 2: rendered as lines
            %   - type >= 2: rendered as scatter plots

            if nargin < 4
                obj.database.addDataset(name, type);
            else
                obj.database.addDataset(name, type, data);
            end
        end

        function obj = addExact(obj, name, exact)
            % ADDEXACT Add exact solution.
            %
            %   obj = addExact(obj, name, exact) registers a new exact
            %   solution with the specified function handle.
            %
            % Inputs:
            %   obj - The Visualizer object
            %   name - Solution name (string or char array)
            %   exact - Exact function handle
            %
            % Outputs:
            %   obj - The Visualizer object

            obj.exacts.(char(name)) = exact;
        end

        function obj = addPlotter(obj, name, strategy)
            % ADDPLOTTER Add tiled plotter.
            %
            %   obj = addPlotter(obj, name, strategy) registers a plotter
            %   with the specified rendering strategy.
            %
            % Inputs:
            %   obj - The Visualizer object
            %   name - Solution name (string)
            %   strategy - Strategy name (string)
            %
            % Outputs:
            %   obj - The Visualizer object

            switch lower(strategy)
                case '1d'
                    nRows = sum(arrayfun(@(x) min(1, x), obj.nComponents));
                    nCols = max(obj.nComponents, [], 1);
                    nAxes = sum(obj.nComponents);
                    nAxesPerComponent = 1;
                    strategy = profilers.visual.Strategy1d();
                case '2d'
                    nRows = sum(obj.nComponents);
                    nCols = length(fieldnames(obj.database.datasets));
                    nAxesPerComponent = nCols;
                    nAxes = nRows * nCols;
                    strategy = profilers.visual.Strategy2d();
                case '3d'
                    nRows = sum(obj.nComponents);
                    nCols = length(fieldnames(obj.database.datasets));
                    nAxesPerComponent = nCols;
                    nAxes = nRows * nCols;
                    strategy = profilers.visual.Strategy3d();
                case 'slice1d'
                    nRows = sum(arrayfun(@(x) min(1, x), obj.nComponents));
                    nCols = max(obj.nComponents, [], 1);
                    nAxes = sum(obj.nComponents);
                    nAxesPerComponent = 1;
                    strategy = profilers.visual.StrategySlice1d();
            end
            obj.plotters.(name) = obj.createPlotter( ...
                nRows, nCols, nAxesPerComponent, nAxes, strategy);
        end

        function obj = initialize(obj)
            % INITIALIZE Set up plotters and figures based on dimensions
            % and components.
            %
            %   obj = initialize(obj) creates the appropriate number of
            %   figures via FigureManager and sets up plotters for the
            %   configured dimensions and components. Also configures
            %   the DataSource with reference geometry for coordinate
            %   generation.
            %
            % Outputs:
            %   obj - The Visualizer object

            %< Reset plotting timeline
            if ~isempty(obj.timeline)
                obj.timeline.reset();
            end

            %< Configure Database with reference geometry
            obj.database.setOrthotope(obj.nDims);
            obj.database.setRefNodes(obj.density);

            if bitand(obj.status, obj.INITIALIZED), return; end

            if obj.isEnabled
                %< Create figures using FigureManager
                plotterNames = fieldnames(obj.plotters);
                nPlotters = length(plotterNames);
                if isempty(obj.figureManager.figureHandles)
                    figureHandles = obj.figureManager.createFigures(nPlotters);
                else
                    figureHandles = obj.figureManager.figureHandles;
                end
    
                %< Attach figures and plotters
                for iPlotter = 1:nPlotters
                    plotterName = plotterNames{iPlotter};
                    plotter = obj.plotters.(plotterName);
                    plotter.setFigure(figureHandles{iPlotter});
                    plotter.draft();
                end
            end

            obj.status = bitor(obj.status, obj.INITIALIZED);
        end

        function obj = plot(obj, state, tDisc)
            % PLOT Generate visualization for current solution state.
            %
            %   obj = plot(obj, state, tDisc) creates plots for the
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

            if nargin < 3, tDisc = []; end

            isSpatialState = isa(state, 'physics.SpatialState');
            core.except.assert(isSpatialState, 'InvalidInput', ...
                'state must be a SpatialState object.');

            space = state.xDisc;
            dofs = state.dofs;

            obj.status = bitor(obj.status, obj.PLOTTING);

            obj.prepare(space);

            if ~isempty(tDisc)
                [V0, n0] = obj.precompute(space, dofs, tDisc);
    
                %< Process each time step
                while obj.timeline.isBefore(tDisc.timeline.next)
                    %< Update all datasets with current solution
                    obj.update(space, dofs, tDisc, V0, n0);
    
                    %< Render all plotters
                    obj.render(obj.timeline.now);
    
                    %< Advance output timeline
                    obj.timeline.advance();
                end
            else
                obj.update(space, dofs);
                obj.render();
            end

            obj.status = bitand(obj.status, bitcmp(obj.PLOTTING));
        end

        function obj = finalize(obj)
            % FINALIZE Clean up all plotters and figures.
            %
            %   obj = finalize(obj) closes all plotters, clears
            %   visualization resources, and closes managed figures.
            %
            % Outputs:
            %   obj - The Visualizer object

            %< Close all plotters
            plotterNames = fieldnames(obj.plotters);
            for i = 1:length(plotterNames)
                plotterName = plotterNames{i};
                obj.plotters.(plotterName).close();
            end
            obj.plotters = struct();

            %< Close all managed figures
            obj.figureManager.closeAll();

            %< Clear database
            obj.database.clear();

            obj.status = bitor(obj.status, obj.FINALIZED);
        end

        function TF = get.isEnabled(obj)
            % GET.ISENABLED Determine if visualization is enabled.
            %
            % Outputs:
            %   TF - True if visualization is enabled (logical)

            TF = false; % Default to disabled
            
            %< Check density initialization and validity
            if isempty(obj.density) || ~isnumeric(obj.density) || any(obj.density <= 0)
                return;
            end
            
            %< Check components structure initialization
            if isempty(obj.components) || ~isstruct(obj.components)
                return;
            else
                names = fieldnames(obj.components);
                if isempty(names)
                    return; % Empty struct
                end
                flag = true;
                for i = 1:length(names)
                    name = names{i};
                    if ~isempty(obj.components.(name))
                        flag = false;
                        break;
                    end
                end
                if flag, return; end
            end

            if bitand(obj.status, obj.INITIALIZED)
                if isempty(obj.plotters) || ~isstruct(obj.plotters)
                    return;
                end
                    
                plotterNames = fieldnames(obj.plotters);
                if isempty(plotterNames)
                    return; % Empty plotters struct
                end
            end
            
            %< All checks passed
            TF = true;
        end

        function TF = get.hasExact(obj)
            % GET.HASEXACT Determine if exact solutions are available.
            %
            % Outputs:
            %   TF - True if exact solutions are available (logical)

            TF = ~isempty(fieldnames(obj.exacts));
        end

        function n = get.nComponents(obj)
            % GET.NCOMPONENTS Get number of components per field.
            %
            % Outputs:
            %   n - Vector of component counts per field

             if isempty(obj.components) || ~isstruct(obj.components)
                n = [];
                return;
            end
            fieldNames = fieldnames(obj.components);
            n = cellfun(@(f) length(obj.components.(f)), fieldNames);
        end

    end

    methods (Access = protected)
        function plotter = createPlotter(obj, nRows, nCols, nAxesPerComponent, nAxes, strategy)
            % CREATEPLOTTER Create configured plotter with strategy.
            %
            % Inputs:
            %   obj - The Visualizer object
            %   nRows - Number of rows in layout
            %   nCols - Number of columns in layout
            %   nAxesPerComponent - Number of axes per component
            %   nAxes - Total number of axes
            %   strategy - Rendering strategy object
            %
            % Outputs:
            %   plotter - Configured TiledPlotter object

            plotter = profilers.visual.TiledPlotter(nRows, nCols);

            %< Bind components to axes
            dofNames = fieldnames(obj.components);
            offset = 0;
            for iField = 1:length(obj.nComponents)
                nFieldComponents = obj.nComponents(iField);
                dofName = dofNames{iField};
                for iComponent = 1:nFieldComponents
                    componentName = sprintf('%s%d', dofName, ...
                        obj.components.(dofName)(iComponent));
                    axisIdx = (offset + iComponent - 1) * nAxesPerComponent + 1;
                    plotter.bind(componentName, axisIdx);
                end
                offset = offset + nFieldComponents;
            end

            %< Configure plotter strategy
            plotter.setStrategy(strategy);

            %< Add axis specifications
            for iAxis = 1:nAxes
                plotter.addAxisSpec(1, 1);
            end
        end

        function obj = prepare(obj, space)
            % PREPARE Set up coordinates and initial datasets.
            %
            % Inputs:
            %   obj - The Visualizer object
            %   space - Spatial discretization object
            %
            % Outputs:
            %   obj - The Visualizer object

            %< Use reference coordinates from Database if available
            if obj.database.hasRefNodes
                phyCoords = space.mesh.collocate(obj.database.refCoords);
                obj.database.setPhysicalCoordinates(phyCoords);
            end
        end

        function [V0, n0] = precompute(obj, space, dofs, tDisc)
            % PRECOMPUTE Compute temporal interpolation data.
            %
            % Inputs:
            %   obj - The Visualizer object
            %   space - Spatial discretization object
            %   dofs - Degrees of freedom structure
            %   tDisc - Temporal discretization object
            %
            % Outputs:
            %   V0 - Cell array of previous solution states
            %   n0 - Number of previous states

            xRef = obj.database.refNodes;
            t0 = tDisc.timeline.now;
            t1 = tDisc.timeline.next;
            dofNames = fieldnames(dofs);

            if t1 > t0
                n0 = min(tDisc.timeline.nSteps, tDisc.timeline.count);
                V0 = cell(1, n0);
                for i = 1:n0
                    W0 = reshape(tDisc.U0{i}, space.nGlobalDofs, []);
                    offset = 0;
                    for j = 1:length(dofNames)
                        dofName = dofNames{j};
                        nComponentsPerDof = size(dofs.(dofName), 2);
                        dof = W0(:, offset+(1:nComponentsPerDof));
                        V0{i}.(dofName) = space.evaluate([], xRef, dof);
                        offset = offset + nComponentsPerDof;
                    end
                end
            else
                V0 = [];
                n0 = 0;
            end
        end

        function obj = update(obj, space, dofs, tDisc, V0, n0)
            % UPDATE Update datasets with current solution.
            %
            % Inputs:
            %   obj - The Visualizer object
            %   space - Spatial discretization object
            %   dofs - Degrees of freedom structure
            %   tDisc - Temporal discretization object
            %   V0 - Previous solution states for interpolation
            %   n0 - Number of previous states
            %
            % Outputs:
            %   obj - The Visualizer object

            if nargin < 4
                tDisc = []; 
                V0 = [];
                n0 = [];
            end

            xRef = obj.database.refNodes;

            if ~isempty(tDisc)
                t = obj.timeline.now;
            end

            dofNames = fieldnames(dofs);

            for i = 1:length(dofNames)
                dofName = dofNames{i};
                if ~isfield(obj.components, dofName)
                    continue;
                end

                component = obj.components.(dofName);
                if isempty(component)
                    continue;
                end

                datasetNames = fieldnames(obj.database.datasets);
                for j = 1:length(datasetNames)
                    datasetName = datasetNames{j};
                    dataset = obj.database.datasets.(datasetName);

                    if dataset.type == dataset.TYPE_FIXED
                        continue; %< Skip fixed datasets
                    end

                    if dataset.type == dataset.TYPE_EXACT
                        %< Exact solution
                        if isfield(obj.exacts, dofName)
                            f = obj.exacts.(dofName);
                            if ~isempty(tDisc)
                                U = space.functionEvaluate([], xRef, f, t);
                            else
                                U = space.functionEvaluate([], xRef, f);
                            end
                            U = U(:, component);
                        else
                            continue; %< No exact solution for this field
                        end
                    else
                        %< Numerical solution
                        U = space.evaluate([], xRef, dofs.(dofName)(:, component));
                        if ~isempty(V0)
                            ht = tDisc.timeline.h;
                            t1 = tDisc.timeline.next;
                            U0 = arrayfun(@(k) V0{k}.(dofName)(:, component), 1:n0, 'Un', 0);
                            U0 = cellfun(@(x) reshape(x, [1, size(x)]), U0, 'Un', 0);
                            n = size(U);
                            U = reshape(U, [1, n]);
                            U = cat(1, U0{:}, U);
                            ti = [t1 - cumsum(ht(n0:-1:1)), t1].';
                            U = interp1(ti, U, t, 'linear', 'extrap');
                            U = reshape(squeeze(U), n);
                        end
                    end

                    %< Store data in dataset using new method
                    obj.database.setDatasetData(datasetName, dofName, U);
                end
            end
        end

        function obj = render(obj, t)
            % RENDER Render all plotters with shared title.
            %
            % Inputs:
            %   obj - The Visualizer object
            %   t - Current time for title generation
            %
            % Outputs:
            %   obj - The Visualizer object

            if nargin < 2, t = []; end

            if ~isempty(t)
                titleText = sprintf('%s at t = %0.2f', obj.titleHead, t);
            else
                titleText = obj.titleHead;
            end
    
            plotterNames = fieldnames(obj.plotters);
            for i = 1:length(plotterNames)
                plotterName = plotterNames{i};
                plotter = obj.plotters.(plotterName);
                plotter.render(titleText, obj.database, []);
            end
        end

    end
end