classdef Visualizer < handle
    %VISUALIZER Own a graphics session; accept already sampled value datasets.
    %   setLayouts accepts Name/Size/Panels and an optional borrowed Figure.
    %   Panel bindings contain Mode, Field, ComponentId and fixed DatasetIds.
    %   render creates graphics lazily. It never samples a numerical method.
    properties (SetAccess = private)
        Layouts
        Style
        FigureHandles
        Axes
        Status = "configured"
        InRender = false
    end
    properties (Dependent)
        IsEnabled
    end
    properties (Access = private)
        States
        Artists
        OwnedFigures
        Visible
    end
    methods
        function obj = Visualizer(layouts, options)
            arguments
                layouts struct = struct([])
                options.style (1, 1) struct = mdemat.tools.visual.plotStyle()
                options.visible (1, 1) string {mustBeMember(options.visible, ["on", "off"])} = "on"
            end
            obj.Layouts = normalizeLayouts(layouts);
            obj.Style = options.style;
            obj.Visible = options.visible;
            obj.clearRegistry();
        end
        function obj = setLayouts(obj, layouts, options)
            arguments
                obj
                layouts struct
                options.style struct = obj.Style
            end
            obj.requireOpen();
            candidate = normalizeLayouts(layouts);
            if ~isequaln(candidate, obj.Layouts)
                report = obj.closeResources(false);
                if ~isempty(report.Vetoed)
                    obj.Status = "finalized";
                    error('mdemat:tools:visual:CloseVetoed', ...
                        'Resolve outstanding owned figures before resetting the session.');
                end
                obj.Layouts = candidate;
                obj.Status = "configured";
            end
            obj.Style = options.style;
        end
        function obj = setStyle(obj, style)
            arguments
                obj
                style (1, 1) struct
            end
            obj.requireOpen();
            if ~all(isfield(style, fieldnames(mdemat.tools.visual.plotStyle())))
                error('mdemat:tools:visual:InvalidStyle', 'Supply a plotStyle value.');
            end
            obj.Style = style;
        end
        function initialize(obj)
            obj.requireOpen();
            previous = get(groot, 'CurrentFigure');
            guard = onCleanup(@() restoreFigure(previous));
            count = numel(obj.Layouts);
            if count == 0, obj.Status = "ready"; return; end
            units = get(groot, 'Units');
            unitsGuard = onCleanup(@() set(groot, 'Units', units));
            set(groot, 'Units', 'pixels');
            screen = get(groot, 'ScreenSize');
            positions = computeFigurePositions(count, screen);
            for k = 1:count
                spec = obj.Layouts(k); name = char(spec.Name);
                valid = isfield(obj.FigureHandles, name) && isgraphics(obj.FigureHandles.(name), 'figure');
                if ~valid
                    if isempty(spec.Figure)
                        fig = figure('Visible', obj.Visible, 'Units', 'pixels', ...
                            'Position', positions(k, :), 'Name', spec.Name);
                        obj.OwnedFigures.(name) = true;
                    elseif isgraphics(spec.Figure, 'figure')
                        fig = spec.Figure;
                        obj.OwnedFigures.(name) = false;
                    else
                        error('mdemat:tools:visual:InvalidBorrowedFigure', 'The borrowed figure was closed.');
                    end
                    obj.FigureHandles.(name) = fig;
                else
                    fig = obj.FigureHandles.(name);
                end
                stateValid = isfield(obj.States, name) && isgraphics(obj.States.(name).Container);
                if stateValid
                    aa = obj.States.(name).Axes;
                    stateValid = all(cellfun(@(a) all(isgraphics(a, 'axes')) && ...
                        all(arrayfun(@(ax) isequal(ancestor(ax, 'figure'), fig) && ...
                        isequal(ancestor(ax, 'uipanel'), obj.States.(name).Container), a)), aa));
                end
                if ~stateValid
                    if isfield(obj.States, name)
                        previousAxes = obj.States.(name).Axes;
                        for a = 1:numel(previousAxes)
                            delete(previousAxes{a}(isgraphics(previousAxes{a})));
                        end
                        if isgraphics(obj.States.(name).Container), delete(obj.States.(name).Container); end
                    end
                    try
                        state = createLayout(fig, spec);
                    catch exception
                        if ~valid && obj.OwnedFigures.(name) && isgraphics(fig), delete(fig); end
                        rethrow(exception);
                    end
                    obj.States.(name) = state;
                    obj.Axes.(name) = state.Axes;
                    obj.Artists.(name) = cell(1, numel(spec.Panels));
                end
            end
            obj.Status = "ready";
        end
        function render(obj, datasets, options)
            arguments
                obj
                datasets struct
                options.time double {mustBeScalarOrEmpty, mustBeFinite} = []
                options.title (1, 1) string = ""
            end
            obj.requireOpen();
            validateDatasets(datasets);
            if obj.InRender
                error('mdemat:tools:visual:ReentrantRender', 'A render is already active.');
            end
            obj.initialize();
            obj.InRender = true;
            previous = get(groot, 'CurrentFigure');
            guard = onCleanup(@() obj.finishRender(previous));
            for k = 1:numel(obj.Layouts)
                spec = obj.Layouts(k); name = char(spec.Name);
                for p = 1:numel(spec.Panels)
                    panel = spec.Panels(p); aa = obj.Axes.(name){p};
                    old = obj.Artists.(name){p};
                    if isempty(panel.RenderFn)
                        delete(old(isgraphics(old)));
                        obj.Artists.(name){p} = gobjects(0, 1);
                        h = mdemat.tools.visual.plotField(aa, datasets, panel.Field, panel.ComponentId, ...
                            mode=panel.Mode, style=obj.Style, datasetIds=panel.DatasetIds, label=panel.Label);
                    else
                        next = arrayfun(@(ax) string(ax.NextPlot), aa);
                        before = arrayfun(@(ax) {ax.Children}, aa);
                        try
                            if strlength(panel.Field) > 0
                                group = mdemat.tools.visual.selectComponent(datasets, panel.Field, ...
                                    panel.ComponentId, datasetIds=panel.DatasetIds);
                            else
                                group = datasets;
                            end
                            h = panel.RenderFn(aa, group);
                            if isstruct(h) && isscalar(h)
                                pieces = struct2cell(h);
                                if ~all(cellfun(@(item) isempty(item) || ...
                                        isobject(item) && all(isgraphics(item(:))), pieces))
                                    error('mdemat:tools:visual:InvalidCallbackArtists', ...
                                        'Handle layers must contain only graphics.');
                                end
                                pieces = cellfun(@(item) item(:), pieces, 'UniformOutput', false);
                                h = vertcat(pieces{:});
                            end
                            if ~obj.ownsArtists(h, aa) || ~obj.declaresNewChildren(h, aa, before)
                                error('mdemat:tools:visual:InvalidCallbackArtists', ...
                                    'Callbacks must return artists in their assigned session axes.');
                            end
                            stale = old(isgraphics(old));
                            stale = stale(~arrayfun(@(item) any(item == h(:)), stale));
                            delete(stale);
                        catch exception
                            delete(aa(isgraphics(aa)));
                            obj.Artists.(name){p} = gobjects(0, 1);
                            rethrow(exception);
                        end
                        for a = 1:numel(aa), if isgraphics(aa(a)), aa(a).NextPlot = char(next(a)); end, end
                    end
                    obj.Artists.(name){p} = h(:);
                end
                text = options.title;
                if ~isempty(options.time), text = text + " at t = " + string(options.time); end
                obj.States.(name).Layout.Title.String = text;
                obj.States.(name).Layout.Title.Interpreter = 'none';
            end
            drawnow;
        end
        function figure = getFigure(obj, name)
            arguments
                obj
                name (1, 1) string
            end
            figure = [];
            if isfield(obj.FigureHandles, name) && isgraphics(obj.FigureHandles.(name), 'figure')
                figure = obj.FigureHandles.(name);
            end
        end
        function report = finalize(obj, options)
            arguments
                obj
                options.force (1, 1) logical = false
            end
            if obj.InRender, error('mdemat:tools:visual:BusySession', 'Wait for rendering before finalizing.'); end
            obj.Status = "finalized";
            report = obj.closeResources(options.force);
        end
        function reset(obj, options)
            arguments
                obj
                options.force (1, 1) logical = false
            end
            if obj.InRender, error('mdemat:tools:visual:BusySession', 'Wait for rendering before resetting.'); end
            report = obj.closeResources(options.force);
            if ~isempty(report.Vetoed)
                error('mdemat:tools:visual:CloseVetoed', 'Owned figures remain open; retry or explicitly force close.');
            end
            obj.Status = "configured";
        end
        function tf = get.IsEnabled(obj)
            tf = ~isempty(obj.Layouts) && obj.Status ~= "finalized";
        end
    end
    methods (Access = private)
        function requireOpen(obj)
            if obj.Status == "finalized"
                error('mdemat:tools:visual:FinalizedSession', 'Call reset before using a finalized session.');
            end
        end
        function finishRender(obj, previous)
            obj.InRender = false;
            restoreFigure(previous);
        end
        function clearRegistry(obj)
            obj.FigureHandles = struct(); obj.Axes = struct(); obj.States = struct();
            obj.Artists = struct(); obj.OwnedFigures = struct();
        end
        function report = closeResources(obj, force)
            report = struct('Vetoed', strings(1, 0));
            failures = {};
            names = fieldnames(obj.FigureHandles);
            for k = 1:numel(names)
                name = names{k}; fig = obj.FigureHandles.(name);
                try
                    if isfield(obj.States, name) && isgraphics(obj.States.(name).Container)
                        delete(obj.States.(name).Container);
                    end
                    if obj.OwnedFigures.(name) && isgraphics(fig, 'figure')
                        if force, delete(fig); else, close(fig); end
                    end
                catch exception
                    failures{end + 1} = exception; %#ok<AGROW>
                end
                if obj.OwnedFigures.(name) && isgraphics(fig, 'figure')
                    report.Vetoed(end + 1) = string(name);
                else
                    obj.FigureHandles = rmfield(obj.FigureHandles, name);
                    obj.OwnedFigures = rmfield(obj.OwnedFigures, name);
                end
            end
            obj.Axes = struct(); obj.States = struct(); obj.Artists = struct();
            if ~isempty(failures)
                exception = MException('mdemat:tools:visual:CloseFailed', ...
                    'Some owned resources could not be closed; their figures remain registered.');
                for j = 1:numel(failures), exception = addCause(exception, failures{j}); end
                throw(exception);
            end
        end
        function tf = declaresNewChildren(~, handles, aa, before)
            tf = true;
            for j = 1:numel(aa)
                children = aa(j).Children;
                for k = 1:numel(children)
                    old = before{j};
                    if ~any(children(k) == old(:)) && ~any(children(k) == handles(:))
                        tf = false; return;
                    end
                end
            end
        end
        function tf = ownsArtists(~, h, aa)
            tf = isempty(h) || isobject(h) && all(isgraphics(h(:)));
            if ~tf, return; end
            for k = 1:numel(h)
                parent = ancestor(h(k), 'axes');
                if isempty(parent) && isappdata(h(k), 'MdeMatOwnerAxes')
                    parent = getappdata(h(k), 'MdeMatOwnerAxes');
                end
                if isempty(parent) || ~any(parent == aa(:)), tf = false; return; end
            end
        end
    end
end

function restoreFigure(previous)
if isempty(previous) || isgraphics(previous, 'figure'), set(groot, 'CurrentFigure', previous); end
end
