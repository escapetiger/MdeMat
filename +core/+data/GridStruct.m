classdef GridStruct < handle
    % GRIDSTRUCT Structure for multidimensional grid data.
    %
    %   GridStruct stores grid coordinates and associated data variables
    %   defined over the grid. Ensures all variables are consistent with
    %   grid dimensions.
    %
    % Examples:
    %   % Create 1D grid
    %   x = linspace(0, 1, 10);
    %   grid = GridStruct(x);
    %   grid.data.temperature = rand(10, 1);
    %
    %   % Create 2D grid
    %   x = linspace(0, 1, 5);
    %   y = linspace(0, 2, 4);
    %   grid = GridStruct({x, y});
    %   grid.data.pressure = rand(20, 1);  % 5*4 = 20 points

    properties
        data   % Structure containing grid data variables
        coords % Grid coordinates
    end

    properties (Dependent)
        nPoints % Total number of grid points
        nDims   % Number of grid dimensions
    end

    methods
        function obj = GridStruct(coords)
            % GRIDSTRUCT Constructor for GridStruct class.
            %
            %   obj = GridStruct() creates empty GridStruct.
            %
            %   obj = GridStruct(coords) creates GridStruct with specified
            %   @a coords and initializes empty data structure.
            %
            % Inputs:
            %   coords - Grid coordinates
            %
            % Outputs:
            %   obj - Constructed GridStruct object

            if nargin < 1, coords = []; end

            obj.data = struct();
            obj.coords = coords;
        end

        function setCoordinates(obj, coords)
            % SETCOORDINATES Set the grid coordinates.
            %
            %   setCoordinates(obj, coords) sets or updates the grid
            %   coordinates @a coords.
            %
            % Inputs:
            %   obj - The GridStruct object
            %   coords - Grid coordinates
            %
            % Outputs:
            %   NULL
            
            if iscell(coords) && numel(coords) == 1
                obj.coords = coords{1};
            else
                obj.coords = coords;
            end
        end

        function n = get.nPoints(obj)
            % GET.NPOINTS Return the total number of grid points.

            if isempty(obj.coords)
                n = 0;
                return;
            end

            if iscell(obj.coords)
                n = prod(cellfun(@numel, obj.coords));
            else
                n = numel(obj.coords);
            end
        end

        function n = get.nDims(obj)
            % GET.NDIMS Return the number of grid dimensions.
            %
            %   n = get.nDims(obj) returns the dimensionality of the grid
            %   based on coordinate specification.
            %
            % Inputs:
            %   obj - The GridStruct object
            %
            % Outputs:
            %   n - Number of dimensions, 0 if no coordinates set

            if isempty(obj.coords)
                n = 0;
            elseif iscell(obj.coords)
                n = length(obj.coords);
            else
                n = 1;
            end
        end
    end
end