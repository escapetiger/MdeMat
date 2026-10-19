classdef Database < handle
    % DATABASE Unified data management for visualization with geometry
    % and grid functionality.
    %
    %   Database provides centralized storage and management of
    %   visualization data using Dataset objects for type-safe data
    %   handling. The class includes reference geometry handling and
    %   coordinate generation, supporting different dataset types
    %   and enabling efficient sharing between multiple plotters with
    %   shared reference coordinate computation.
    %
    %   The class manages both reference coordinates (for evaluation on
    %   canonical domains) and physical coordinates (actual spatial
    %   locations), enabling consistent coordinate systems across all
    %   visualization components.
    %
    % Examples:
    %   % Create database with geometry
    %   database = Database();
    %   database.setOrthotope(2);
    %   database.setRefNodes(15); % 15×15 evaluation grid
    %   
    %   % Add datasets and set physical coordinates
    %   database.addDataset('numerical', 2);
    %   database.addDataset('exact', 1);
    %   
    %   % Parse data for visualization
    %   plotterData = database.parse('u1');
    %
    % See Also:
    %   profilers.visual.Dataset, profilers.visual.Visualizer,
    %   profilers.visual.TiledPlotter

    properties
        datasets % Structure containing named Dataset objects
        geometry % Reference geometry element defining the domain
        refNodes % Reference evaluation nodes for rendering (nDims × nPhyPoints)
        refCoords % Cell array of reference coordinate vectors for each dimension
    end

    properties (Dependent)
        nDims % Number of spatial dimensions (scalar)
        hasGeometry % Flag indicating if geometry is configured (logical)
        hasRefNodes % Flag indicating if reference nodes are configured (logical)
    end

    methods
        function obj = Database()
            % DATABASE Constructor for Database.
            %
            %   obj = Database() creates an empty database ready
            %   for configuration and data storage.
            %
            % Outputs:
            %   obj - Constructed Database object

            obj.datasets = struct();
            obj.geometry = [];
            obj.refNodes = [];
            obj.refCoords = {};
        end

        function obj = setOrthotope(obj, nDims)
            % SETORTHOTOPE Configure unit orthotope geometry for coordinate generation.
            %
            %   obj = setOrthotope(obj, nDims) creates and sets a unit
            %   orthotope geometry in the specified number of dimensions
            %   for generating evaluation nodes and coordinates.
            %
            % Inputs:
            %   obj - The Database object
            %   nDims - Number of spatial dimensions (positive integer)
            %
            % Outputs:
            %   obj - The Database object

            obj.geometry = core.geometry.Orthotope.unit(nDims);
        end

        function obj = setRefNodes(obj, n)
            % SETREFNODES Configure evaluation nodes for rendering.
            %
            %   obj = setRefNodes(obj, n) sets up a uniform grid of evaluation
            %   nodes within the reference geometry bounds. These nodes are
            %   used for function evaluation and data visualization. Requires
            %   geometry to be set first via setOrthotope().
            %
            % Inputs:
            %   obj - The Database object
            %   n - Number of nodes per dimension (scalar or vector)
            %
            % Outputs:
            %   obj - The Database object

            core.except.assert(~isempty(obj.geometry), 'NoGeometry', ...
                'Geometry must be set before configuring nodes via setOrthotope()');

            d = obj.geometry.nDims;
            a = obj.geometry.lower;
            b = obj.geometry.upper;
            
            %< Handle scalar input for uniform grids
            if isscalar(n)
                n = repmat(n, 1, d);
            end
            
            %< Create coordinate vectors for each dimension
            h = (b - a) ./ (n + 1);
            obj.refCoords = arrayfun(@(i) (a(i)+h(i)):h(i):(b(i)-h(i)), 1:d, 'Un', false);

            %< Generate tensor product grid
            [nodeGrids{1:d}] = ndgrid(obj.refCoords{:});
            obj.refNodes = reshape(cat(d+1, nodeGrids{:}), [], d).';
        end

        function obj = setPhysicalCoordinates(obj, phyCoords)
            % SETPHYSICALCOORDINATES Configure physical coordinates.
            %
            %   obj = setPhysicalCoordinates(obj, n) sets up physical
            %   coordinates for all datasets that are not of fixed type.
            %
            % Inputs:
            %   obj - The Database object
            %   phyCoords - Physical coordinates (vector)
            %
            % Outputs:
            %   obj - The Database object

            datasetNames = fieldnames(obj.datasets);
            for iDataset = 1:length(datasetNames)
                datasetName = datasetNames{iDataset};
                dataset = obj.datasets.(datasetName);
                if dataset.type ~= dataset.TYPE_FIXED
                    dataset.setData('x', phyCoords);
                end
            end
        end

        function obj = addDataset(obj, name, type, data)
            % ADDDATASET Register a new dataset with specified type.
            %
            %   obj = addDataset(obj, name, type) creates a new Dataset
            %   object with the specified name and type.
            %
            %   obj = addDataset(obj, name, type, data) creates a new
            %   Dataset object with the specified name, type and initial
            %   data structure.
            %
            % Inputs:
            %   obj - The Database object
            %   name - Unique name for the dataset (string)
            %   type - Dataset type (0=fixed, 1=exact, 2=numerical)
            %   data - Initial data structure (optional)
            %
            % Outputs:
            %   obj - The Database object

            name = char(name);
            
            %< Create Dataset object
            if nargin < 4 || isempty(data)
                dataset = profilers.visual.Dataset(type);
            else
                dataset = profilers.visual.Dataset(type, data);
            end
            obj.datasets.(name) = dataset;
        end

        function varDataset = parse(obj, varName)
            % PARSE Parse data for rendering from datasets.
            %
            %   varDataset = parse(obj, varName) extracts and formats data
            %   for the specified variable from all datasets, organizing
            %   coordinates, values, and data types for visualization.
            %   Always uses physical coordinates for consistency across
            %   visualization components.
            %
            % Inputs:
            %   obj - The Database object
            %   varName - Name of the variable to parse (e.g., 'u1', 'pressure2')
            %
            % Outputs:
            %   varDataset - Structure mapping dataset names to parsed data

            %< Parse variable name to extract field and component index
            tokens = regexp(varName, '^([a-zA-Z_]+)(\d+)', 'tokens');
            core.except.assert(~isempty(tokens), 'InvalidInput', ...
                'Variable name does not match expected format');
            
            fieldName = tokens{1}{1};
            varIdx = str2double(tokens{1}{2});

            varDataset = struct();
            datasetNames = fieldnames(obj.datasets);
            
            for i = 1:length(datasetNames)
                datasetName = datasetNames{i};
                dataset = obj.datasets.(datasetName);

                %< Skip if field not present in this dataset
                if ~dataset.hasField(fieldName)
                    continue;
                end

                fieldData = dataset.getData(fieldName);
                
                %< Extract the specific component
                if size(fieldData, 2) < varIdx
                    continue;
                end
                componentData = fieldData(:, varIdx);

                %< Always use physical coordinates for rendering
                phyCoords = dataset.getData('x');
                if iscell(phyCoords)
                    coords = phyCoords(1:obj.nDims);
                else
                    coords = cell(1, obj.nDims);
                    for iDim = 1:obj.nDims
                        coords{iDim} = unique(phyCoords(iDim, :));
                    end
                end

                %< Reshape data if needed
                if obj.nDims == 1
                    n = [1, length(coords{1})];
                else
                    n = cellfun(@length, coords);
                end
                
                if length(componentData) == prod(n)
                    reshapedData = reshape(componentData, n(:).');
                else
                    reshapedData = componentData;
                end

                %< Create data structure with type information
                dataStruct = struct();
                dataStruct.u = reshapedData;
                dataStruct.type = dataset.type;
                dataStruct.nDims = obj.nDims;
                
                %< Add coordinate information
                for iDim = 1:obj.nDims
                    coordName = sprintf('x%d', iDim);
                    dataStruct.(coordName) = coords{iDim};
                end

                varDataset.(datasetName) = dataStruct;
            end
        end

        function obj = clear(obj)
            % CLEAR Remove all datasets and reset state.
            %
            %   obj = clear(obj) removes all datasets and resets
            %   coordinate and geometry state.
            %
            % Inputs:
            %   obj - The Database object
            %
            % Outputs:
            %   obj - The Database object

            obj.datasets = struct();
            obj.geometry = [];
            obj.refNodes = [];
            obj.refCoords = {};
        end

        function dataset = getDataset(obj, name)
            % GETDATASET Get Dataset object by name.
            %
            %   dataset = getDataset(obj, name) retrieves the Dataset
            %   object with the specified name.
            %
            % Inputs:
            %   obj - The Database object
            %   name - Name of the dataset (string or char)
            %
            % Outputs:
            %   dataset - Dataset object (empty if not found)

            name = char(name);
            if isfield(obj.datasets, name)
                dataset = obj.datasets.(name);
            else
                dataset = [];
            end
        end

        function obj = setDatasetData(obj, datasetName, fieldName, fieldData)
            % SETDATASETDATA Set data for a field in a specific dataset.
            %
            %   obj = setDatasetData(obj, datasetName, fieldName, fieldData)
            %   sets the data for the specified field in the named dataset.
            %
            % Inputs:
            %   obj - The Database object
            %   datasetName - Name of the dataset (string or char)
            %   fieldName - Name of the field (string or char)
            %   fieldData - Data array for the field
            %
            % Outputs:
            %   obj - The Database object

            dataset = obj.getDataset(datasetName);
            if ~isempty(dataset)
                dataset.setData(fieldName, fieldData);
            else
                core.except.verify(false, 'DatasetNotFound', ...
                    'Dataset ''%s'' not found', char(datasetName));
            end
        end

        function fieldData = getDatasetData(obj, datasetName, fieldName)
            % GETDATASETDATA Get data for a field from a specific dataset.
            %
            %   fieldData = getDatasetData(obj, datasetName, fieldName)
            %   retrieves the data for the specified field from the named dataset.
            %
            % Inputs:
            %   obj - The Database object
            %   datasetName - Name of the dataset (string or char)
            %   fieldName - Name of the field (string or char)
            %
            % Outputs:
            %   fieldData - Data array for the field (empty if not found)

            dataset = obj.getDataset(datasetName);
            if ~isempty(dataset)
                fieldData = dataset.getData(fieldName);
            else
                fieldData = [];
            end
        end

        function n = get.nDims(obj)
            % GET.NDIMS Get number of spatial dimensions.
            %
            %   n = get.nDims(obj) returns the number of spatial
            %   dimensions based on available coordinates or geometry.
            %
            % Outputs:
            %   n - Number of dimensions (scalar)

            if ~isempty(obj.refCoords)
                n = length(obj.refCoords);
            elseif ~isempty(obj.geometry)
                n = obj.geometry.nDims;
            else
                n = 0;
            end
        end

        function TF = get.hasGeometry(obj)
            % GET.HASGEOMETRY Check if geometry is configured.
            %
            % Outputs:
            %   TF - True if geometry is set (logical)

            TF = ~isempty(obj.geometry);
        end

        function TF = get.hasRefNodes(obj)
            % GET.HASREFNODES Check if reference nodes are configured.
            %
            % Outputs:
            %   TF - True if reference nodes are available (logical)

            TF = ~isempty(obj.refNodes) && ~isempty(obj.refCoords);
        end
    end
end