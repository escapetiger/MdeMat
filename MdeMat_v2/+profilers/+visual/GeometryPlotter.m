classdef GeometryPlotter < profilers.visual.TiledPlotter
    % GEOMETRYPLOTTER Specialized plotter for reference geometry
    % visualization.
    %
    %   GeometryPlotter extends TiledPlotter to provide specialized
    %   functionality for visualizing data on reference geometry elements.
    %   This class automatically configures rendering nodes and coordinates
    %   based on geometry bounds, simplifying visualization workflows that
    %   don't require full finite element space setup.
    %
    %   The plotter is particularly useful for visualizing data on
    %   canonical domains (unit squares, cubes, etc.) and for prototyping
    %   visualization strategies before applying them to complex finite
    %   element meshes.
    %
    % Examples:
    %   % Create plotter for unit square with 2×3 layout
    %   geometry = core.geometry.Orthotope.unit(2);
    %   plotter = GeometryPlotter(geometry, 2, 3);
    %   
    %   % Configure rendering and strategy
    %   plotter.setNodes(15); % 15×15 evaluation grid
    %   plotter.setStrategy(profilers.visual.Strategy2d('contour'));
    %   
    %   % Set up layout and render
    %   plotter.addAxisSpec(1, 1); % Add axis specifications
    %   plotter.draft('TileSpacing', 'compact');
    %   plotter.render('Function Visualization', datasets, styles);
    %
    % See also:
    %   profilers.visual.TiledPlotter, core.geometry.Orthotope,
    %   profilers.visual.Strategy

    properties (SetAccess = protected, GetAccess = public)
        geometry % Reference geometry element defining the domain
        nodes    % Reference evaluation nodes for rendering (nDims×nPoints)
        coords   % Cell array of coordinate vectors for each dimension
    end

    methods
        function obj = GeometryPlotter(geometry, nRows, nCols)
            % GEOMETRYPLOTTER Constructor for GeometryPlotter.
            %
            %   obj = GeometryPlotter(geometry, nRows, nCols) creates a
            %   geometry plotter for the specified reference geometry
            %   with a tiled layout of the given dimensions.
            %
            % Inputs:
            %   geometry - Reference geometry object (e.g., Orthotope)
            %   nRows - Number of rows in the tiled layout (positive integer)
            %   nCols - Number of columns in the tiled layout (positive integer)
            %
            % Outputs:
            %   obj - Constructed GeometryPlotter object

            obj@profilers.visual.TiledPlotter(nRows, nCols);
            obj.geometry = geometry;
            obj.nodes = {};
            obj.coords = {};
        end

        function obj = setNodes(obj, n)
            % SETNODES Configure evaluation nodes for rendering.
            %
            %   obj = setNodes(obj, n) sets up a uniform grid of evaluation
            %   nodes within the reference geometry bounds. These nodes are
            %   used for function evaluation and data visualization.
            %
            % Inputs:
            %   obj - The GeometryPlotter object
            %   n - Number of nodes per dimension (scalar or vector)
            %
            % Outputs:
            %   obj - The GeometryPlotter object

            d = obj.geometry.nDims;
            a = obj.geometry.lower;
            b = obj.geometry.upper;
            
            %< Handle scalar input for uniform grids
            if isscalar(n)
                n = repmat(n, 1, d);
            end
            
            %< Create coordinate vectors for each dimension
            h = (b - a) ./ (n + 1);
            obj.coords = arrayfun(@(i) (a(i)+h(i)):h(i):(b(i)-h(i)), 1:d, 'Un', false);

            %< Generate tensor product grid
            [obj.nodes{1:d}] = ndgrid(obj.coords{:});
            obj.nodes = reshape(cat(d+1, obj.nodes{:}), [], d).';
        end
    end
end