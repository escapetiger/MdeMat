classdef Compiler < handle
    % COMPILER Base class for all symbolic compilers.
    %
    %   Compiler provides a unified interface for generating, compiling,
    %   and storing basis matrices from various symbolic function types. It
    %   manages a centralized repository of pre-computed basis matrices for
    %   efficient reuse across the core package.
    %
    %   The compiler creates function handle matrices for basis functions
    %   and their derivatives, stores them with metadata for later retrieval,
    %   and provides utilities for managing the compiled function repository.
    %
    % See Also:
    %   core.symbolic.InterpolationPolynomialBasisCompiler,
    %   core.symbolic.OrthogonalPolynomialBasisCompiler,
    %   core.symbolic.FourierBasisCompiler

    properties (Constant)
        METADATA_DIR = 'config/functions' % Directory for storing compiled matrices
        FILE_EXTENSION = '.mat' % File extension for compiled matrices
    end

    properties
        variables % Cell array of symbolic variables
        funcs % Object array of SymbolicFunction instances
    end

    methods
        function obj = Compiler(variables)
            % COMPILER Constructor for Compiler.
            %
            %   obj = Compiler() creates a new Compiler object with the
            %   default variable {sym('x')}.
            %
            %   obj = Compiler(variables) creates a new Compiler object with
            %   the specified symbolic variables and ensures the metadata
            %   directory exists for storing compiled basis matrices.
            %
            % Inputs:
            %   variables - (Optional) Cell array of symbolic variables, default: {sym('x')}
            %
            % Outputs:
            %   obj - Constructed Compiler object

            if ~exist(obj.METADATA_DIR, 'dir')
                [success, message] = mkdir(obj.METADATA_DIR);
                core.except.assert(success, 'DirectoryCreationFailed', ...
                    'Failed to create metadata directory %s: %s', ...
                    obj.METADATA_DIR, message);
            end

            if nargin < 1
                variables = {sym('x')};
            end
            obj.variables = variables;
            obj.funcs = [];
        end

        function filename = compile(obj, fileName, basisName, maxOrder)
            % COMPILE Compiles polynomial basis matrix with derivatives.
            %
            %   filename = compile(obj, fileName, basisName, maxOrder) creates
            %   function handle matrices for basis functions and their
            %   derivatives up to the specified order. Saves the result
            %   with metadata for later retrieval.
            %
            % Inputs:
            %   obj - The Compiler object
            %   fileName - Filename for the compiled matrix (without extension)
            %   basisName - Name/type of the basis functions for documentation
            %   maxOrder - Maximum derivative order to compute (non-negative integer)
            %
            % Outputs:
            %   filename - Generated filename for the stored matrix

            core.except.assert(~isempty(obj.funcs), ...
                'InvalidState', 'No basis functions set. Call setFunctions first.');

            core.except.assert(maxOrder >= 0 && isscalar(maxOrder) && ...
                maxOrder == floor(maxOrder), 'InvalidInput', ...
                'Maximum order must be a non-negative integer.');

            handles = obj.createHandles(maxOrder);

            metadata = obj.createMetadata(lower(basisName));
            metadata.basisName = basisName;
            metadata.numFunctions = length(obj.funcs);
            metadata.maxOrder = maxOrder;
            metadata.size = [length(obj.funcs), maxOrder + 1];
            metadata.description = sprintf('%s basis matrix with %d funcs, max order %d', ...
                basisName, length(obj.funcs), maxOrder);
            metadata.handles = handles;

            obj.save(fileName, metadata);

            fprintf('Compiled %s basis matrix: %s\n', basisName, fileName);
            filename = fileName;
        end

        function metadata = load(obj, filename)
            % LOAD Loads a previously compiled basis matrix.
            %
            %   metadata = load(obj, filename) retrieves a compiled basis
            %   matrix from the metadata directory and returns the metadata
            %   structure containing function handles and compilation information.
            %
            % Inputs:
            %   obj - The Compiler object
            %   filename - Filename of the stored matrix (without extension)
            %
            % Outputs:
            %   metadata - Metadata structure containing handles and compilation info

            fullPath = obj.getFullPath(filename);

            core.except.assert(exist(fullPath, 'file') == 2, ...
                'FileNotFound', 'Basis matrix file not found: %s', fullPath);

            data = load(fullPath);
            metadata = data.metadata;
        end

        function fileList = list(obj)
            % LIST Lists all compiled basis matrices.
            %
            %   fileList = list(obj) returns a list of all available
            %   compiled basis matrices in the metadata directory.
            %
            % Inputs:
            %   obj - The Compiler object
            %
            % Outputs:
            %   fileList - Cell array of compiled matrix filenames (without extension)

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

        function obj = delete(obj, filename)
            % DELETE Deletes a compiled basis matrix.
            %
            %   obj = delete(obj, filename) removes the specified compiled
            %   basis matrix file from the metadata directory.
            %
            % Inputs:
            %   obj - The Compiler object
            %   filename - File name without extension
            %
            % Outputs:
            %   obj - The Compiler object

            fullPath = obj.getFullPath(filename);

            if exist(fullPath, 'file')
                try
                    delete(fullPath);
                    fprintf('Deleted compiled matrix: %s\n', filename);
                catch ME
                    core.except.verify(0, 'DeletionFailed', ...
                        'Failed to delete file %s: %s', fullPath, ME.message);
                end
            else
                core.except.verify(0, 'FileNotFound', ...
                    'File not found: %s', fullPath);
            end
        end

        function obj = clear(obj)
            % CLEAR Deletes all compiled basis matrices.
            %
            %   obj = clear(obj) removes all compiled basis matrix files
            %   from the metadata directory. Useful for cleanup or when
            %   starting fresh with a new compilation session.
            %
            % Inputs:
            %   obj - The Compiler object
            %
            % Outputs:
            %   obj - The Compiler object

            fileList = obj.list();

            for i = 1:length(fileList)
                obj.delete(fileList{i});
            end

            fprintf('Cleared %d compiled matrices.\n', length(fileList));
        end
    end

    methods (Abstract)
        % SETFUNCTIONS Abstract method for setting basis functions.
        obj = setFunctions(obj, varargin)
    end

    methods (Access = protected)
        function obj = save(obj, filename, metadata)
            % SAVE Saves metadata to file.
            %
            %   obj = save(obj, filename, metadata) stores the compiled
            %   basis matrix metadata to a MAT file in the metadata
            %   directory using MATLAB v7.3 format for large data
            %   compatibility.
            %
            % Inputs:
            %   obj - The Compiler object
            %   filename - Filename for the saved matrix
            %   metadata - Metadata structure to save
            %
            % Outputs:
            %   obj - The Compiler object

            fullPath = obj.getFullPath(filename);

            try
                save(fullPath, 'metadata', '-v7.3');
            catch ME
                core.except.assert(0, 'SaveFailed', ...
                    'Failed to save basis matrix to %s: %s', fullPath, ME.message);
            end
        end

        function fullPath = getFullPath(obj, filename)
            % GETFULLPATH Constructs full file path from filename.
            %
            %   fullPath = getFullPath(obj, filename) creates the complete
            %   file path by combining the metadata directory with the
            %   filename and appropriate extension if not already present.
            %
            % Inputs:
            %   obj - The Compiler object
            %   filename - Base filename
            %
            % Outputs:
            %   fullPath - Full file path with extension

            if ~endsWith(filename, obj.FILE_EXTENSION)
                filename = [filename, obj.FILE_EXTENSION];
            end

            fullPath = fullfile(obj.METADATA_DIR, filename);
        end

        function B = createHandles(obj, maxOrder)
            % CREATEHANDLES Creates basis matrix from stored functions.
            %
            %   B = createHandles(obj, maxOrder) converts symbolic basis
            %   functions to function handles and computes derivatives up
            %   to the specified order. Returns a cell matrix of function
            %   handles for efficient numerical evaluation.
            %
            % Inputs:
            %   obj - The Compiler object
            %   maxOrder - Maximum derivative order
            %
            % Outputs:
            %   B - Cell matrix of function handles [nFunctions x (maxOrder+1)]

            core.except.assert(~isempty(obj.funcs), ...
                'InvalidState', 'No basis functions available.');

            n = length(obj.funcs);
            m = maxOrder + 1;
            x = obj.variables;

            B = cell(n, m);
            S = cell(n, m);

            for i = 1:n
                S{i, 1} = obj.funcs(i).expression;
            end

            for i = 1:n
                expr = S{i, 1};
                for j = 1:(m - 1)
                    S{i, j+1} = simplify(diff(expr, x, j));
                end
            end

            for i = 1:n
                for j = 1:m
                    expr = S{i, j};
                    if ~isempty(expr)
                        B{i, j} = matlabFunction(expr, 'Vars', x);
                    end
                end
            end
        end

        function metadata = createMetadata(obj, basisType)
            % CREATEMETADATA Creates base metadata structure.
            %
            %   metadata = createMetadata(obj, basisType) generates the
            %   base metadata structure with common fields for all compiled
            %   basis matrices including timestamp, version information,
            %   and compilation details.
            %
            % Inputs:
            %   obj - The Compiler object
            %   basisType - Type of basis function
            %
            % Outputs:
            %   metadata - Base metadata structure

            metadata = struct( ...
                'type', basisType, ...
                'variables', {obj.variables}, ...
                'timestamp', datetime('now'), ...
                'matlabVersion', version(), ...
                'generator', class(obj) ...
                );
        end
    end
end