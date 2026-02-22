classdef UniformGrid < approx.mesh.Grid
    % UNIFORMGRID Uniform multidimensional structured grid.
    %
    %   UniformGrid represents a structured grid with uniform spacing
    %   defined by the number of elements per dimension and domain bounds.
    %   This class implements memory-efficient tensor-based storage with
    %   constant element sizes throughout the domain.
    %
    %   Uniform grids have constant element spacing within each dimension.
    %   The bounding box format is [a1, b1, a2, b2, ..., ad, bd] where
    %   [ai, bi] defines the interval for dimension i.
    %
    % See also:
    %   approx.mesh.Grid, approx.mesh.NonuniformGrid,
    
    methods
        function obj = UniformGrid(bbox, resolution)
            % UNIFORMGRID Create uniform grid from bounding box and resolution.
            %
            %   obj = UniformGrid(bbox, resolution) creates a uniform grid
            %   with specified bounding box and element count per dimension.
            
            arguments
                bbox{mustBeNumeric}
                resolution{mustBePositive, mustBeInteger}
            end
            
            bbox = bbox(:).';
            nx = resolution(:).';
            nd = length(bbox) / 2;
            a = bbox(1:2:end);
            b = bbox(2:2:end);
            vertices = arrayfun(@(i) linspace(a(i), b(i), nx(i)+1), 1:nd, 'Un', 0);
            obj@approx.mesh.Grid(vertices);
        end
        
        function newObj = refine(obj, nLevels)
            % REFINE Create refined uniform grid preserving spacing
            % patterns.
            %
            %   newObj = refine(obj, nLevels) creates a refined uniform
            %   grid by uniformly subdividing each element while preserving
            %   the relative spacing characteristics of the original grid.
            
            arguments
                obj approx.mesh.UniformGrid
                nLevels{mustBeNonnegative, mustBeInteger} = 0
            end
            
            nx = obj.Resolution;
            k = 2^nLevels;
            a = cellfun(@(x) min(x), obj.Vertices);
            b = cellfun(@(x) max(x), obj.Vertices);
            bbox = [a(:).'; b(:).'];
            newObj = approx.mesh.UniformGrid(bbox, nx.*k);
        end
    end
end
