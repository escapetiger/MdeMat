classdef Dataset < handle
    % DATASET Data container for visualization datasets.
    %
    %   Dataset provides a simple container for visualization data with
    %   type information. This class encapsulates the data structure and
    %   type classification used throughout the visualization system.
    %
    %   The @a Type property determines how the data is rendered:
    %   - Type 0 (FIXED): Fixed or reference data
    %   - Type 1 (EXACT): Exact analytical solution data
    %   - Type 2 (NUMERICAL): Numerical solution data
    %
    % See also:
    %   physics.visual.Database
    
    properties
        Type % Type of dataset (0: FIXED, 1: EXACT, 2: NUMERIC)
        Data % Struct containing field data
    end
    
    methods
        function obj = Dataset(options)
            % DATASET Construct an instance of Dataset.
            %
            %   obj = Dataset() creates an empty dataset with default type and
            %   empty data structure.
            %
            %   obj = Dataset(type=type) creates a dataset with specified @a type
            %   and empty data structure.
            %
            %   obj = Dataset(type=type, data=data) creates a dataset with
            %   specified @a type and @a data structure.
            
            arguments
                options.type {mustBeMember(options.type, {'fixed', 'exact', 'numeric'})} = 'numeric'
                options.data struct = struct()
            end
            
            obj.Type = options.type;
            obj.Data = options.data;
        end
        
        function obj = setData(obj, fieldName, fieldData)
            % SETDATA Set data for a specific field.
            %
            %   obj = setData(obj, fieldName, fieldData) sets the @a fieldData
            %   for the specified @a fieldName in the dataset.
            
            arguments
                obj physics.visual.Dataset
                fieldName {mustBeTextScalar}
                fieldData {mustBeNumeric}
            end
            
            obj.Data.(char(fieldName)) = fieldData;
        end
        
        function fieldData = getData(obj, fieldName)
            % GETDATA Get data for a specific field.
            %
            %   fieldData = getData(obj, fieldName) retrieves the @a fieldData
            %   for the specified @a fieldName from the dataset. Returns empty
            %   array if field does not exist.
            
            arguments
                obj physics.visual.Dataset
                fieldName {mustBeTextScalar}
            end
            
            fieldName = char(fieldName);
            if isfield(obj.Data, fieldName)
                fieldData = obj.Data.(fieldName);
            else
                fieldData = [];
            end
        end
        
        function TF = hasField(obj, fieldName)
            % HASFIELD Check if dataset contains a specific field.
            %
            %   TF = hasField(obj, fieldName) returns true if the
            %   dataset contains data for the specified @a fieldName.
            
            arguments
                obj physics.visual.Dataset
                fieldName {mustBeTextScalar}
            end
            
            TF = isfield(obj.Data, char(fieldName));
        end
        
        function fieldNames = getFieldNames(obj)
            % GETFIELDNAMES Get list of available field names.
            %
            %   fieldNames = getFieldNames(obj) returns a cell array
            %   of field names available in this dataset.
            
            arguments
                obj physics.visual.Dataset
            end
            
            fieldNames = fieldnames(obj.Data);
        end
        
        function obj = clearData(obj)
            % CLEARDATA Remove all data from the dataset.
            %
            %   obj = clearData(obj) removes all field data from
            %   the dataset while preserving the @a Type.
            
            arguments
                obj physics.visual.Dataset
            end
            
            obj.Data = struct();
        end
    end
end