classdef Dataset < handle
    % DATASET Data container for visualization datasets.
    %
    %   Dataset provides a simple container for visualization data with
    %   type information. This class encapsulates the data structure and
    %   type classification used throughout the visualization system.
    %
    %   The type property determines how the data is rendered:
    %   - Type 0 (FIXED): Fixed or reference data
    %   - Type 1 (EXACT): Exact analytical solution data
    %   - Type 2 (NUMERICAL): Numerical solution data
    %
    % Examples:
    %   % Create dataset with numerical data
    %   dataset = Dataset(2, struct('u', numericalSolution));
    %   
    %   % Create dataset with exact solution
    %   exactData = struct('u', exactSolution, 'v', exactVelocity);
    %   dataset = Dataset(1, exactData);
    %   
    %   % Access properties
    %   dataType = dataset.type;
    %   fieldData = dataset.data.u;
    %
    % See Also:
    %   profilers.visual.DataSource
    
    properties (Constant)
        % Dataset type constants
        TYPE_FIXED = 0; % Fixed/reference data
        TYPE_EXACT = 1; % Exact solution data
        TYPE_NUMERICAL = 2; % Numerical solution data
    end

    properties
        type % Dataset type (0=fixed, 1=exact, 2=numerical)
        data % Data structure containing field arrays
        phyCoords % Physical coordinate array for all datasets (cell array)
    end

    methods
        function obj = Dataset(type, data)
            % DATASET Constructor for Dataset.
            %
            %   obj = Dataset() creates an empty dataset with default
            %   type and empty data structure.
            %
            %   obj = Dataset(type) creates a dataset with specified
            %   type and empty data structure.
            %
            %   obj = Dataset(type, data) creates a dataset with
            %   specified type and data structure.
            %
            % Inputs:
            %   type - Dataset type (scalar, default: 2)
            %   data - Data structure (struct, default: empty struct)
            %
            % Outputs:
            %   obj - Constructed Dataset object

            if nargin < 1
                obj.type = 2; %< Default to numerical type
            else
                obj.type = type;
            end

            if nargin < 2
                obj.data = struct();
            else
                obj.data = data;
            end

            obj.phyCoords = [];
        end

        function obj = setData(obj, fieldName, fieldData)
            % SETDATA Set data for a specific field.
            %
            %   obj = setData(obj, fieldName, fieldData) sets the data
            %   for the specified field name in the dataset.
            %
            % Inputs:
            %   obj - The Dataset object
            %   fieldName - Name of the field (string or char)
            %   fieldData - Data array for the field
            %
            % Outputs:
            %   obj - The Dataset object

            obj.data.(char(fieldName)) = fieldData;
        end

        function fieldData = getData(obj, fieldName)
            % GETDATA Get data for a specific field.
            %
            %   fieldData = getData(obj, fieldName) retrieves the data
            %   for the specified field name from the dataset.
            %
            % Inputs:
            %   obj - The Dataset object
            %   fieldName - Name of the field (string or char)
            %
            % Outputs:
            %   fieldData - Data array for the field (empty if not found)

            fieldName = char(fieldName);
            if isfield(obj.data, fieldName)
                fieldData = obj.data.(fieldName);
            else
                fieldData = [];
            end
        end

        function TF = hasField(obj, fieldName)
            % HASFIELD Check if dataset contains a specific field.
            %
            %   TF = hasField(obj, fieldName) returns true if the
            %   dataset contains data for the specified field name.
            %
            % Inputs:
            %   obj - The Dataset object
            %   fieldName - Name of the field (string or char)
            %
            % Outputs:
            %   TF - True if field exists (logical)

            TF = isfield(obj.data, char(fieldName));
        end

        function fieldNames = getFieldNames(obj)
            % GETFIELDNAMES Get list of available field names.
            %
            %   fieldNames = getFieldNames(obj) returns a cell array
            %   of field names available in this dataset.
            %
            % Inputs:
            %   obj - The Dataset object
            %
            % Outputs:
            %   fieldNames - Cell array of field name strings

            fieldNames = fieldnames(obj.data);
        end

        function obj = clearData(obj)
            % CLEARDATA Remove all data from the dataset.
            %
            %   obj = clearData(obj) removes all field data from
            %   the dataset while preserving the type.
            %
            % Inputs:
            %   obj - The Dataset object
            %
            % Outputs:
            %   obj - The Dataset object

            obj.data = struct();
        end
    end
end