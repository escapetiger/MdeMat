classdef SymbolicFunction < handle
    % SYMBOLICFUNCTION Base class for all symbolic functions.
    %
    %   SymbolicFunction provides a foundation for creating and
    %   manipulating symbolic mathematical expressions with a consistent
    %   interface. This class serves as a base for implementing specialized
    %   symbolic function types while maintaining uniform evaluation and
    %   conversion capabilities.
    %
    % Examples:
    %   % Create a symbolic function with default variable
    %   obj = SymbolicFunction();
    %
    %   % Create with custom variables
    %   vars = {sym('x'), sym('y')};
    %   obj = SymbolicFunction(vars);
    %
    % See also:
    %   sym, syms, matlabFunction
    
    properties (Access = public)
        expression  % The symbolic expression
        variables   % Cell array of symbolic variables
    end
    
    methods
        function obj = SymbolicFunction(variables)
            % SYMBOLICFUNCTION Constructor for SymbolicFunction.
            %
            %   SymbolicFunction() creates a new SymbolicFunction object
            %   with the default variable sym('x').
            %
            %   SymbolicFunction(variables) creates a new SymbolicFunction
            %   object with the specified symbolic variables.
            %
            % Inputs:
            %   variables - Symbol(s) (optional, default: {sym('x')})
            %
            % Outputs:
            %   obj - Constructed SymbolicFunction object
            
            if nargin < 1, variables = {sym('x')}; end

            if ~iscell(variables), variables = {variables}; end

            obj.variables = variables;
            obj.expression = [];
        end
        
        function result = evaluate(obj, values)
            % EVALUATE Evaluates the symbolic expression at specified
            % values.
            %
            %   result = evaluate(obj, values) converts the symbolic
            %   expression to a function handle and evaluates it at the
            %   given numeric values. The number of values must match the
            %   number of variables in the function.
            %
            % Inputs:
            %   obj - The SymbolicFunction object
            %   values - Numeric value(s) (scalar, vector or cell) 
            %
            % Outputs:
            %   result - Numeric result of the evaluation
            
            if ~iscell(values) && numel(obj.variables) > 1
                core.except.assert( ...
                    numel(values) == numel(obj.variables), ...
                    'DimensionMismatch', ...
                    'Number of values must match number of variables.');
                values = mat2cell(values);
            end
            
            try
                f = matlabFunction(obj.expression, 'Vars', obj.variables);
                if iscell(values)
                    result = f(values{:});
                else
                    result = f(values);
                end
            catch ME
                core.except.assert(0, 'EvaluationFailed', ...
                    'Failed to evaluate expression: %s', ME.message);
            end
        end
        
        function result = toString(obj)
            % TOSTRING Converts the symbolic expression to a string.
            %
            %   result = toString(obj) returns a character string
            %   representation of the symbolic expression using MATLAB's
            %   char function.
            %
            % Inputs:
            %   obj - The SymbolicFunction object
            %
            % Outputs:
            %   result - Character string representation of the expression
            
            try
                result = char(obj.expression);
            catch ME
                core.except.assert(0, 'ConversionFailed', ...
                    'Failed to convert expression to string: %s', ME.message);
            end
        end
        
        function result = toFunctionHandle(obj)
            % TOFUNCTIONHANDLE Converts the symbolic expression to a
            % function handle.
            %
            %   result = toFunctionHandle(obj) creates a MATLAB function
            %   handle from the symbolic expression that can be used for
            %   efficient numerical evaluation.
            %
            % Inputs:
            %   obj - The SymbolicFunction object
            %
            % Outputs:
            %   result - MATLAB function handle
            
            try
                result = matlabFunction(obj.expression, 'Vars', obj.variables);
            catch ME
                core.except.assert(0, 'ConversionFailed', ...
                    'Failed to convert expression to function handle: %s', ME.message);
            end
        end
        
        function result = toLatex(obj)
            % TOLATEX Converts the symbolic expression to a LaTeX string.
            %
            %   result = toLatex(obj) generates a LaTeX representation of
            %   the symbolic expression suitable for mathematical
            %   typesetting and documentation.
            %
            % Inputs:
            %   obj - The SymbolicFunction object
            %
            % Outputs:
            %   result - LaTeX string representation of the expression

            try
                result = latex(obj.expression);
            catch ME
                core.except.assert(0, 'ConversionFailed', ...
                    'Failed to convert expression to LaTeX string: %s', ...
                    ME.message);
            end
        end
    end
end