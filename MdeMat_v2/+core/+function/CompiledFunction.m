classdef CompiledFunction < core.function.Function
    % COMPILEDFUNCTION Base class for all compiled functions.
    %
    %   CompiledFunction provides a bridge between symbolic function
    %   definitions and efficient numerical evaluation by automatically
    %   loading precompiled function handles and metadata from the
    %   configuration directory.
    %
    %   The class supports automatic metadata discovery, caching, and
    %   efficient function evaluation using precompiled function handles
    %   stored in MAT files.
    %
    % See Also:
    %   core.function.Function, core.function.DifferentiableFunction

    properties (Constant, Access = protected)
        METADATA_DIR = 'config/functions' % Directory containing metadata files
        FILE_EXTENSION = '.mat' % File extension for metadata files
    end

    properties
        metadata % Loaded metadata structure containing function handles
    end

    properties (Dependent)
        hasMetadata % True if metadata is loaded and available
        metasource % Cell array of available metadata filenames
    end

    methods
        function obj = CompiledFunction(nDims, nCodims)
            % COMPILEDFUNCTION Constructor for the CompiledFunction class.
            %
            %   obj = CompiledFunction() creates a CompiledFunction with
            %   default dimensions (both nDims and nCodims are 0).
            %   
            %   obj = CompiledFunction(nDims) creates a CompiledFunction
            %   with specified input dimension and default output dimension
            %   (0).
            %   
            %   obj = CompiledFunction(nDims, nCodims) creates a
            %   CompiledFunction with specified input and output
            %   dimensions.
            %
            % Inputs:
            %   nDims - Number of input dimensions (optional, default:0)
            %   nCodims - Number of output dimensions (optional, default: 0)
            %
            % Outputs:
            %   obj - The created CompiledFunction object

            obj@core.function.Function(nDims, nCodims);
            obj.metadata = struct();
        end

        function TF = get.hasMetadata(obj)
            % GET.HASMETADATA Getter for hasMetadata dependent property.
            %
            %   TF = get.hasMetadata(obj) returns true if the object has
            %   loaded metadata with valid function handles.
            %
            % Inputs:
            %   obj - The CompiledFunction object
            %
            % Outputs:
            %   TF - True if metadata is loaded and contains handles

            TF = ~isempty(fieldnames(obj.metadata)) && ...
                isfield(obj.metadata, 'handles') && ...
                ~isempty(obj.metadata.handles);
        end

        function fileList = get.metasource(obj)
            % GET.METASOURCE Getter for metasource dependent property.
            %
            %   fileList = get.metasource(obj) returns a cell array of
            %   available metadata filenames in the metadata directory.
            %
            % Inputs:
            %   obj - The CompiledFunction object
            %
            % Outputs:
            %   fileList - Cell array of available metadata filenames

            if ~exist(obj.METADATA_DIR, 'dir')
                fileList = {};
                return;
            end

            files = dir(fullfile(obj.METADATA_DIR, ['*', obj.FILE_EXTENSION]));
            fileList = cell(length(files), 1);

            for i = 1:length(files)
                [~, name, ~] = fileparts(files(i).name);
                fileList{i} = name;
            end
        end

        function obj = load(obj, filename)
            % LOAD Load compiled metadata from specified file.
            %
            %   obj = load(obj, filename) loads precompiled function
            %   metadata from the specified file and updates the object's
            %   function dimensions based on the loaded metadata.
            %
            % Inputs:
            %   obj - The CompiledFunction object
            %   filename - Metadata filename
            %
            % Outputs:
            %   obj - The CompiledFunction object

            core.except.assert( ...
                ischar(filename) || isstring(filename), ...
                'InvalidInput', ...
                'Filename must be a string or char array.');

            filename = char(filename);
            fullPath = obj.getFullPath(filename);

            core.except.assert(exist(fullPath, 'file') == 2, ...
                'FileNotFound', ...
                'Metadata file not found: %s', fullPath);

            data = load(fullPath, 'metadata');

            core.except.assert(isfield(data, 'metadata'), ...
                'InvalidFile', ...
                'Metadata file must contain a metadata field.');

            obj.metadata = data.metadata;
        end

        function obj = autoLoad(obj)
            % AUTOLOAD Automatically find and load appropriate metadata.
            %
            %   obj = autoLoad(obj) automatically determines the
            %   appropriate metadata file for this function type and loads
            %   it. The filename is determined by the abstract autoFilename
            %   method.
            %
            % Inputs:
            %   obj - The CompiledFunction object
            %
            % Outputs:
            %   obj - The CompiledFunction object

            filename = obj.autoFilename();

            core.except.assert(~isempty(filename) && ...
                (ischar(filename) || isstring(filename)), ...
                'MetadataNotFound', ...
                'No appropriate metadata file found for %s', class(obj));

            obj = obj.load(filename);
        end
    end

    methods (Access = protected)
        function Y = evaluateImpl(obj, X)
            % EVALUATEIMPL Implementation of function evaluation.
            %
            %   Y = evaluateImpl(obj, X) evaluates the compiled function
            %   at the specified points using the loaded metadata handles.
            %
            % Inputs:
            %   obj - The CompiledFunction object
            %   X - Coordinates (nDims x nPoints)
            %
            % Outputs:
            %   Y - Function values (nCodims x nPoints)

            core.except.assert(obj.hasMetadata, ...
                'NoMetadata', 'No compiled metadata is loaded.');

            handles = obj.metadata.handles;
            Y = zeros(obj.nCodims, size(X, 2));

            Z = mat2cell(X, ones(1, size(X, 1)), size(X, 2));

            for j = 1:size(handles, 1)
                if ~isempty(handles{j, 1})
                    Y(j, :) = handles{j, 1}(Z{:});
                end
            end
        end
        
        function dY = derivativeImpl(obj, X, r)
            % DERIVATIVEIMPL Implementation of function derivative
            % evaluation.
            %
            %   dY = derivativeImpl(obj, X, r) computes derivatives of the
            %   compiled function at the specified points using precompiled
            %   derivative handles from the loaded metadata.
            %
            % Inputs:
            %   obj - The DifferentiableCompiledFunction object
            %   X - Coordinates (nDims x nPoints)
            %   r - Derivative order vector (1 x nDims)
            %
            % Outputs:
            %   dY - Derivative values (nCodims x nPoints)

            core.except.assert(obj.hasMetadata, ...
                'NoMetadata', 'No compiled metadata is loaded.');

            handles = obj.metadata.handles;
            
            maxOrder = size(handles, 2) - 1;
            
            derivativeOrder = sum(r);
            
            core.except.assert(derivativeOrder <= maxOrder, ...
                'InvalidInput', ...
                'Derivative order %d exceeds maximum available order %d.', ...
                derivativeOrder, maxOrder);

            columnIndex = derivativeOrder + 1;
            dY = zeros(obj.nCodims, size(X, 2));
            
            for j = 1:size(handles, 1)
                if ~isempty(handles{j, columnIndex})
                    dY(j, :) = handles{j, columnIndex}(X);
                end
            end
        end
  
        function fullPath = getFullPath(obj, filename)
            % GETFULLPATH Construct full file path from filename.
            %
            %   fullPath = getFullPath(obj, filename) constructs the full
            %   path to a metadata file, adding the file extension if
            %   needed.
            %
            % Inputs:
            %   obj - The CompiledFunction object
            %   filename - Base filename string
            %
            % Outputs:
            %   fullPath - Full path to the metadata file

            if ~endsWith(filename, obj.FILE_EXTENSION)
                filename = [filename, obj.FILE_EXTENSION];
            end

            fullPath = fullfile(obj.METADATA_DIR, filename);
        end
    end

    methods (Abstract, Access = protected)
        % AUTOFILENAME Abstract method to generate appropriate metadata
        % filename.
        filename = autoFilename(obj)
    end
end