classdef Database < handle
    % DATABASE Simplified data container for visualization datasets.
    %
    %   Database provides centralized storage and management of
    %   visualization data using Dataset objects. The class focuses
    %   purely on data management and provides grouping functionality
    %   for visualization components.
    %
    % See also:
    %   physics.visual.Dataset, physics.visual.Visualizer,
    %   physics.visual.TiledPlotter
    
    properties
        Datasets % Structure containing named Dataset objects
    end
    
    properties
        NDatasets % Number of datasets in the database
        DatasetNames % Cell array of dataset names
        IsEmpty % True if the database has no datasets
    end
    
    methods
        function obj = Database()
            % DATABASE Construct an instance of Database.
            %
            %   obj = Database() creates an empty database ready for
            %   dataset storage.
            
            obj.Datasets = struct();
        end
        
        function n = get.NDatasets(obj)
            % GET.NDATASETS Get the number of datasets in the database.
            
            n = length(fieldnames(obj.Datasets));
        end

        function names = get.DatasetNames(obj)
            % GET.DATASETNAMES Get the names of all datasets in the database.
            
            names = fieldnames(obj.Datasets);
        end

        function tf = get.IsEmpty(obj)
            % GET.ISEMPTY Check if the database is empty.
            
            tf = obj.NDatasets == 0;
        end
        
        function obj = setDataset(obj, name, dataset)
            % SETDATASET Store a dataset with the specified name.
            %
            %   obj = setDataset(obj, name, dataset) stores the @a dataset
            %   with the specified @a name in the database.
            
            arguments
                obj physics.visual.Database
                name {mustBeTextScalar}
                dataset physics.visual.Dataset
            end
            
            obj.Datasets.(char(name)) = dataset;
        end
        
        function dataset = getDataset(obj, name)
            % GETDATASET Retrieve a dataset by name.
            %
            %   dataset = getDataset(obj, name) retrieves the Dataset
            %   object with the specified @a name.
            
            arguments
                obj physics.visual.Database
                name {mustBeTextScalar}
            end
            
            name = char(name);
            if isfield(obj.Datasets, name)
                dataset = obj.Datasets.(name);
            else
                dataset = [];
            end
        end
        
        function database = groupBy(obj, variable)
            % GROUPBY Create a new database with processed field data.
            %
            %   database = groupBy(obj, variable) extracts and formats data
            %   for the specified @a variable from all datasets, organizing
            %   coordinates, values, and data types for visualization.
            %   Returns a new Database instance with processed datasets.
            
            arguments
                obj physics.visual.Database
                variable {mustBeTextScalar}
            end
            
            tokens = regexp(variable, '^([a-zA-Z_]+)(\d+)', 'tokens');
            core.except.assert(~isempty(tokens), 'InvalidInput', ...
                'Variable name does not match expected format');
            
            fieldName = tokens{1}{1};
            varIdx = str2double(tokens{1}{2});
            
            database = physics.visual.Database();
            datasetNames = fieldnames(obj.Datasets);
            
            for i = 1:length(datasetNames)
                datasetName = datasetNames{i};
                dataset = obj.Datasets.(datasetName);
                
                if ~dataset.hasField(fieldName)
                    continue;
                end
                
                fieldData = dataset.getData(fieldName);
                
                if size(fieldData, 2) < varIdx
                    continue;
                end
                componentData = fieldData(:, varIdx);
                
                x = dataset.getData('x');
                if ismatrix(x)
                    nd = size(x, 1);
                    x = arrayfun(@(d) unique(x(d, :)), 1:nd, 'Un', 0);
                else
                    nd = length(x);
                end
                
                if nd == 1
                    n = [1, length(x{1})];
                else
                    n = cellfun(@length, x);
                end
                
                if length(componentData) == prod(n)
                    u = reshape(componentData, n(:).');
                else
                    u = componentData;
                end
                
                data = struct();
                data.x = x;
                data.u = u;
                
                newDataset = physics.visual.Dataset(Type=dataset.Type, Data=data);
                database.setDataset(datasetName, newDataset);
            end
        end
    end
end