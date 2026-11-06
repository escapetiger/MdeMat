function h = plot(mesh, options)
%PLOT Draw supported mesh geometry without storing a graphics session in it.
%   Triangulation inherits this method. Filled elements are simplices; general
%   polytopes require an explicit supported topology and are never convexified.
arguments
    mesh (1, 1) mdemat.data.mesh.PolytopalMesh
    options.axes = []
    options.vertexColor = 'b'
    options.edgeColor = 'k'
    options.faceColor = 'none'
    options.elementColor = 'none'
    options.boundaryColor = 'r'
    options.vertexSize (1, 1) double {mustBePositive, mustBeFinite} = 36
    options.edgeWidth (1, 1) double {mustBePositive, mustBeFinite} = 1
    options.alpha (1, 1) double {mustBeNonnegative, mustBeLessThanOrEqual(options.alpha, 1), mustBeFinite} = 0.3
    options.showVertices (1, 1) logical = true
    options.showEdges (1, 1) logical = true
    options.showFaces (1, 1) logical = false
    options.showElements (1, 1) logical = false
    options.showBoundary (1, 1) logical = true
    options.useFastRendering (1, 1) logical = true
    options.replace logical {mustBeScalarOrEmpty} = []
    options.maxEdges (1, 1) double {mustBeNonnegative} = Inf
end
h = struct('Vertices', gobjects(0, 1), 'Edges', gobjects(0, 1), ...
    'Faces', gobjects(0, 1), 'Elements', gobjects(0, 1), ...
    'Boundary', gobjects(0, 1), 'Legend', gobjects(0, 1));
v = mesh.Vertices; nd = mesh.NDims;
if isempty(v), return; end
if ~isnumeric(v) || ~isreal(v) || ~ismatrix(v) || any(~isfinite(v(:)))
    error('mdemat:data:mesh:PolytopalMesh:InvalidVertices', 'Vertices must be finite real rows.');
end
if nd <= 3 && mesh.NVerticesPerElement ~= nd + 1
    error('mdemat:data:mesh:PolytopalMesh:UnsupportedTopology', 'Only explicit simplex elements are supported.');
end
edges = zeros(0, 2);
if options.showEdges
    graph = mesh.graphify(); edges = graph.Edges;
    if size(edges, 1) > options.maxEdges
        error('mdemat:data:mesh:PolytopalMesh:EdgeBudgetExceeded', 'Mesh edges exceed maxEdges.');
    end
end
explicit = ~isempty(options.axes);
if explicit, ax = options.axes; else, ax = gca; end
if ~isscalar(ax) || ~isgraphics(ax, 'axes')
    error('mdemat:data:mesh:PolytopalMesh:InvalidAxes', 'Supply one valid axes.');
end
next = ax.NextPlot; previous = get(groot, 'CurrentFigure');
guard = onCleanup(@() restoreAxes(ax, next, previous));
try
    replace = options.replace;
    if isempty(replace), replace = ~explicit && ~ishold(ax); end
    if replace, cla(ax); end
    ax.NextPlot = 'add';
    if nd > 3
        centered = v - mean(v, 1);
        [~, ~, directions] = svd(centered, 'econ');
        v = centered * directions(:, 1:min(3, size(directions, 2)));
        if size(v, 2) < 3, v(:, end + 1:3) = 0; end
        warning('mdemat:data:mesh:PolytopalMesh:DimensionReduction', ...
            'Visualizing %dD mesh by a three-dimensional PCA projection.', nd);
    end
    base = mdemat.data.mesh.internal.plotGraph(v, edges, axes=ax, ...
        vertexColor=options.vertexColor, edgeColor=options.edgeColor, ...
        vertexSize=options.vertexSize, edgeWidth=options.edgeWidth, ...
        showVertices=options.showVertices, showEdges=options.showEdges, showReference=false);
    h.Vertices = base.Vertices; h.Edges = base.Edges;
    if nd <= 3
        if options.showElements && ~strcmp(options.elementColor, 'none') && mesh.NElements > 0
            switch nd
                case 1
                    h.Elements = mdemat.data.mesh.internal.drawEdges(ax, v, mesh.Elements, options.elementColor, 3);
                case 2
                    h.Elements = patch(ax, 'Faces', mesh.Elements, 'Vertices', v, ...
                        'FaceColor', options.elementColor, 'FaceAlpha', options.alpha, 'EdgeColor', 'none');
                case 3
                    if options.useFastRendering
                        faces = mesh.Faces(mesh.Boundary, :);
                    else
                        local = nchoosek(1:4, 3);
                        faces = reshape(permute(reshape(mesh.Elements(:, local.'), ...
                            mesh.NElements, 3, 4), [1, 3, 2]), [], 3);
                    end
                    if ~isempty(faces)
                        h.Elements = patch(ax, 'Faces', faces, 'Vertices', v, ...
                            'FaceColor', options.elementColor, 'FaceAlpha', options.alpha, 'EdgeColor', 'none');
                    end
            end
        end
        if options.showFaces && ~strcmp(options.faceColor, 'none') && mesh.NFaces > 0
            switch nd
                case 1
                    indices = unique(mesh.Faces(:));
                    h.Faces = scatter(ax, v(indices, 1), zeros(numel(indices), 1), ...
                        options.vertexSize, options.faceColor, 'filled');
                case 2
                    h.Faces = mdemat.data.mesh.internal.drawEdges(ax, v, mesh.Faces, options.faceColor, options.edgeWidth);
                case 3
                    h.Faces = patch(ax, 'Faces', mesh.Faces, 'Vertices', v, ...
                        'FaceColor', options.faceColor, 'FaceAlpha', options.alpha, 'EdgeColor', 'none');
            end
        end
        if options.showBoundary && ~strcmp(options.boundaryColor, 'none') && mesh.NBoundaryFaces > 0
            faces = mesh.Faces(mesh.Boundary, :);
            switch nd
                case 1
                    indices = unique(faces(:));
                    h.Boundary = scatter(ax, v(indices, 1), zeros(numel(indices), 1), ...
                        100, options.boundaryColor, 'filled', 'MarkerEdgeColor', 'k');
                case 2
                    h.Boundary = mdemat.data.mesh.internal.drawEdges(ax, v, faces, options.boundaryColor, 1.5 * options.edgeWidth);
                case 3
                    h.Boundary = patch(ax, 'Faces', faces, 'Vertices', v, ...
                        'FaceColor', 'none', 'EdgeColor', options.boundaryColor, 'LineWidth', 1.5 * options.edgeWidth);
            end
        end
    end
    labels = fieldnames(h); entries = gobjects(0, 1); names = strings(1, 0);
    for k = 1:numel(labels)
        handles = h.(labels{k});
        if ~isempty(handles), entries(end + 1, 1) = handles(1); names(end + 1) = string(labels{k}); end
    end
    if ~isempty(entries) && isempty(ax.Legend)
        h.Legend = legend(ax, entries, names, 'Interpreter', 'none');
        setappdata(h.Legend, 'MdeMatOwnerAxes', ax);
    end
    axis(ax, 'equal'); grid(ax, 'on');
    if nd >= 3, view(ax, 3); else, view(ax, 2); end
    title(ax, sprintf('%dD PolytopalMesh: %d vertices, %d elements, %d faces (%d boundary)', ...
        nd, mesh.NVertices, mesh.NElements, mesh.NFaces, mesh.NBoundaryFaces));
catch exception
    parts = struct2cell(h);
    for k = 1:numel(parts), delete(parts{k}(isgraphics(parts{k}))); end
    rethrow(exception);
end
end

function restoreAxes(ax, next, previous)
if isgraphics(ax), ax.NextPlot = next; end
if isempty(previous) || isgraphics(previous, 'figure'), set(groot, 'CurrentFigure', previous); end
end
