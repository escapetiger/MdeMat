classdef UniformGrid < approx.mesh.Grid
    % UNIFORMGRID Uniform multidimensional structured grid.
    %
    %   UniformGrid represents a structured grid with uniform spacing
    %   defined by the number of elements per dimension and domain bounds.
    %   This class implements memory-efficient tensor-based storage with
    %   constant element sizes throughout the domain.
    %
    % Examples:
    %   % 2D uniform grid: 10×10 elements on unit square
    %   grid2D = UniformGrid([10, 10], [0, 0, 1, 1]);
    %   
    %   % 3D uniform grid: 5×8×12 elements on custom domain
    %   grid3D = UniformGrid([5, 8, 12], [-1, 0, -2, 1, 2, 3]);
    %   
    %   % Query grid properties
    %   h = grid2D.measure;  % Minimum element size
    %   nElems = grid2D.nTotalElements;
    %
    % Notes:
    %   Uniform grids have constant element spacing within each dimension.
    %   The bounding box format is [a1, b1, a2, b2, ..., ad, bd] where
    %   [ai, bi] defines the interval for dimension i.
    %
    % See also:
    %   approx.mesh.Grid, approx.mesh.NonuniformGrid,
    
    methods
        function obj = UniformGrid(bbox, resolution)
            % UNIFORMGRID Create uniform grid from bounding box and
            % resolution.
            %
            %   obj = UniformGrid(bbox, resolution) creates a uniform grid
            %   with specified bounding box and element count per
            %   dimension.
            %
            % Inputs:
            %   bbox - Bounding box coordinates (1 × 2*nDims vector)
            %   resolution - Element count per dimension (1 × nDims vector)
            %
            % Outputs:
            %   obj - Constructed UniformGrid object

            bbox = bbox(:).';
            n = resolution(:).';
            d = length(bbox) / 2;
            a = bbox(1:2:end);
            b = bbox(2:2:end);
            vertices = arrayfun(@(i) linspace(a(i), b(i), n(i)+1), 1:d, 'Un', 0);
            obj@approx.mesh.Grid(vertices);
        end

        function newObj = refine(obj, nLevels)
            % REFINE Create refined uniform grid preserving spacing
            % patterns.
            %
            %   newObj = refine(obj, nLevels) creates a refined unfirom
            %   grid by uniformly subdividing each element while preserving
            %   the relative spacing characteristics of the original grid.
            %
            % Inputs:
            %   obj - The UniformGrid object
            %   nLevels - Number of refinement levels (non-negative integer)
            %
            % Outputs:
            %   newObj - Refined UniformGrid object

            n = obj.resolution;
            k = 2^nLevels;
            a = cellfun(@(x) min(x), obj.vertices);
            b = cellfun(@(x) max(x), obj.vertices);
            bbox = [a(:).'; b(:).'];
            newObj = approx.mesh.UniformGrid(bbox, n .* k);
        end
    end
end

