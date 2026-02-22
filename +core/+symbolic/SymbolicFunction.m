classdef SymbolicFunction < handle
    % SYMBOLICFUNCTION Base class for all symbolic functions.
    %
    %   SymbolicFunction represents a symbolic function with expressions
    %   and variables. It provides the foundation for all symbolic basis
    %   functions with compilation capabilities and dimension properties.
    %
    % See also:
    %   core.symbolic.LagrangeBasis, core.symbolic.LegendreBasis,
    %   core.symbolic.ChebyshevBasis
    
    properties (Constant)
        MetaExt = '.mat' % Metadata file extension
    end
    
    properties
        Expressions % Object array of symbolic expressions
        Variables % Object array of symbolic variables
        Handles % Cell array of function handles
    end
    
    properties (Dependent)
        NDims % Number of input dimensions
        NCodims % Number of output dimensions
    end
    
    methods
        function obj = SymbolicFunction(variables)
            % SYMBOLICFUNCTION Construct an instance of SymbolicFunction.
            %
            %   obj = SymbolicFunction(variables) creates a symbolic
            %   function with the specified @a variables.
            
            arguments
                variables {mustBeA(variables, 'sym')} = []
            end
            
            obj.Variables = variables(:);
            obj.Expressions = sym.empty(0, 1);
            obj.Handles = {};
        end
        
        function obj = compile(obj, maxDerivOrder)
            % COMPILE Compile symbolic function as function handles.
            %
            %   obj = compile(obj) creates function handles with default
            %   derivative order 0.
            %
            %   obj = compile(obj, maxOrder) creates function handles up to
            %   the specified @a maxDerivOrder.
            
            arguments
                obj core.symbolic.SymbolicFunction
                maxDerivOrder {mustBeNonnegative, mustBeInteger} = 0
            end
            
            isCompilable = obj.NDims == 1 || maxDerivOrder == 0;
            core.except.assert(isCompilable, 'NotCompilable', ...
                'Derivatives for multivariate functions not supported.');
            
            obj.Handles = cell(obj.NCodims, maxDerivOrder+1);
            
            for i = 1:obj.NCodims
                expr = obj.Expressions(i);
                
                for j = 0:maxDerivOrder
                    if j == 0
                        obj.Handles{i, j+1} = matlabFunction( ...
                            expr, 'Vars', num2cell(obj.Variables));
                        continue;
                    end
                    
                    derivExpr = diff(expr, obj.Variables(1), j);
                    obj.Handles{i, j+1} = matlabFunction( ...
                        derivExpr, 'Vars', num2cell(obj.Variables));
                end
            end
        end
        
        function obj = save(obj, fileName)
            % SAVE Save symbolic function metadata to file.
            %
            %   obj = save(fileName) saves the symbolic function metadata
            %   to the specified @a fileName.
            
            arguments
                obj core.symbolic.SymbolicFunction
                fileName {mustBeTextScalar}
            end
            
            [~, ~, ext] = fileparts(fileName);
            if isempty(ext)
                fileName = [fileName, obj.MetaExt];
            end
            
            metadata = struct();
            metadata.nDims = obj.NDims;
            metadata.nCodims = obj.NCodims;
            metadata.variables = obj.Variables;
            metadata.expressions = obj.Expressions;
            metadata.handles = obj.Handles;
            metadata.compileTime = datetime('now');
            
            save(fileName, 'metadata');
        end
        
        function obj = load(obj, fileName)
            % LOAD Load symbolic function metadata from file and update
            % object.
            %
            %   obj = load(fileName) loads the symbolic function metadata
            %   from the specified @a fileName and updates the current
            %   object properties with the loaded data.
            
            arguments
                obj core.symbolic.SymbolicFunction
                fileName {mustBeTextScalar}
            end
            
            [~, ~, ext] = fileparts(fileName);
            if isempty(ext)
                fileName = [fileName, obj.MetaExt];
            end
            
            core.except.assert(isfile(fileName), 'FileNotFound', ...
                'File not found: %s', fileName);
            
            loadedData = load(fileName);
            core.except.assert(isfield(loadedData, 'metadata'), ...
                'InvalidFile', 'Invalid metadata file: %s', fileName);
            
            metadata = loadedData.metadata;
            
            obj.Variables = metadata.variables;
            obj.Expressions = metadata.expressions;
            if isfield(metadata, 'handles')
                obj.Handles = metadata.handles;
            else
                obj.Handles = {};
            end
        end
    end
    
    methods
        function n = get.NDims(obj)
            % GET.NDIMS Returns the number of input dimensions.
            n = length(obj.Variables);
        end
        
        function n = get.NCodims(obj)
            % GET.NCODIMS Returns the number of output dimensions.
            n = length(obj.Expressions);
        end
    end
end