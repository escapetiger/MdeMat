%% Graphics: value samples and one session owner
% Graph data owns its plot interface; it never owns a figure or axes:
%
%   graph = mdemat.data.mesh.Graph([0,0;1,1],[1,2]);
%   fig = figure; ax = axes(fig);
%   h = plot(graph,axes=ax);
%
% h.Vertices and h.Edges are graphics handles. PolytopalMesh (including
% Triangulation) and SeparableGraph provide the same method convention.
% An explicit axes is additive unless replace=true. Omitting axes preserves
% ordinary interactive hold-off replacement. maxVertices preserves factor
% endpoints only when the budget permits it; otherwise it raises an error.
%
%% A field is a value, not a live Dataset object
% Values has components-by-points layout. ComponentIds names original rows:
%
%   [X,Y] = ndgrid([-1,0,2],[-2,1,4]);
%   field = struct('Values',[X(:).';Y(:).'],'ComponentIds',[2,4]);
%   field.Axes = {[-1,0,2],[-2,1,4]};
%   samples.RUN = struct('Type',"numeric",'Data',struct('u',field));
%   selected = mdemat.tools.visual.selectComponent(samples,"u",4);
%
% Data.u.Values(2,:) is component 4. Neither selectComponent nor the graphics
% layer infers row 4. Axes follows ndgrid column-major order. Alternatively,
% Coordinates is dimension-by-points, paired with Values without sorting.
% A point view is never silently converted to a tensor or interpolated.
% Physical fields named x, v or Values are safe: geometry is inside the view.
% Type (fixed/exact/numeric) provides default styling, not a sampling policy.
%
%% Draw once, or own a session
%   style = mdemat.tools.visual.plotStyle(symmetricColorLimits=true);
%   h = mdemat.tools.visual.plotField(ax,samples,"u",4,mode="2d",style=style);
%   native = namedargs2cell(mdemat.tools.visual.curveStyle(style,2,1));
%   % plot(ax,x,y,native{:}) uses the same color and pattern selection.
%
% plotField returns newly owned artists to its caller. It preserves existing
% unrelated artists and temporary NextPlot/current-figure state. A 2-D tensor
% contains sampling nodes; surfaces retain nonuniform positions and all rows.
% Missing values may be NaN. Infinite values are rejected; complex display
% requires an explicit real/imag/abs/phase choice. Only the selected display
% row is densified when sparse; other components and State are not expanded. Three-dimensional volume
% and diagonal interpolation require two nodes per tensor axis.
%
%   panel = struct('Mode',"2d",'Field',"u",'ComponentId',4,...
%       'DatasetIds',"RUN");
%   layout = struct('Name',"field",'Size',[1,1],'Panels',panel);
%   session = mdemat.tools.visual.Visualizer(layout,style=style);
%   session.render(samples,time=0,title="Example");
%   fig = session.getFigure("field");
%   session.finalize();
%
% A layout has Name, Size=[rows columns], Panels, and optional borrowed Figure.
% A panel has Mode, Field, ComponentId and fixed DatasetIds. Optional Tiles
% are row-major tile starts; Spans is one [rows columns] row per axis. Omitted
% starts are packed before any graphics are allocated. TileSpacing and Padding
% retain the native tiledlayout appearance options. 1d and slice1d overlay
% fixed data sets on one axis; 2d/3d use one axis per fixed dataset ID.
% Missing samples leave that binding empty and do not shift axes or colors.
%
%% Resource lifetime and callbacks
% Visualizer is the only graphics owner. It does not retain numerical State,
% reference functions, spatial density, a timeline or samples. Submit sampled
% values explicitly. FigureHandles, Axes and Layouts are read-only queries.
% Copying a numeric sample struct does not create an implicitly updated view.
%
% For custom panels, set RenderFn=@(axes,group)... and Mode="custom". A callback
% must return its new/updated artists from the assigned session-owned axes,
% as a graphics vector or a scalar struct of handle layers returned by data.plot.
% Field may be empty to receive all samples. Callback failure discards its
% exclusive axes; another render recreates them. Returning an unrelated artist
% is rejected. Side effects outside the assigned axes are the callback's own.
%
% finalize is terminal; render requires explicit reset afterwards. Close
% callbacks are respected: report=session.finalize() names vetoed figures,
% which remain registered. reset refuses while they remain open; resolve the
% callback or explicitly use reset(force=true). Borrowed figures are never
% closed. A user-closed owned window is recreated on the next render.
%
% The former Dataset, Database, Strategy, FigureManager, TiledPlotter and mesh
% Plotter classes have been replaced. There is no compatibility object tree.
