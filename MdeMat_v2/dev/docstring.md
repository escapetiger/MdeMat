# MdeMat Docstring Template

#### Function Documentation Template

```matlab
function [output1, output2] = functionName(input1, input2, varargin)
% FUNCTIONNAME One-line summary of function purpose.
%
%   [output1, output2] = functionName(input1, input2) computes something. 
%   Detailed description of syntax and action of the function.
%
% Inputs:
%   input1 - Description of first input argument
%   input2 - Description of second input argument
%
% Outputs:
%   output1 - Description of first output
%   output2 - Description of second output
%   varargin - Additional arguments
%     'Parameter' - 
%
% Examples:
%   % Basic usage
%   result = functionName(input);
%
%   % Advanced usage with parameters
%   [out1, out2] = functionName(in1, in2, 'Parameter', value);
%
% Notes: (optional)
%   Any special considerations, limitations, or implementation details.
%
% See Also:
% 	package.module.func1, package.module.func2
end
```

#### Class Documentation Template

```matlab
classdef ClassName < ParentClass
    % CLASSNAME One-line summary of class purpose.
    %
    %   Extended description of class purpose and functionality.
    %   This section should outline the major capabilities and use cases.
    %
    % Examples:
    %   obj = ClassName(arg1, arg2);
    %   result = obj.method1(input);
    %
    % Notes: (optional)
    %   Any special considerations, limitations, or implementation details.
    %
    % See also:
    %   package.module.cls1, package.module.func1
    
    properties
        prop1 % Description of the first property
        
        prop2 % Description of the second property
    end
    
    methods
        function obj = ClassName(arg1, arg2)
            % CLASSNAME Constructor for ClassName.
            %
            %   obj = ClassName(obj, arg1) does ... 
            %
            %   obj = ClassName(obj, arg1, arg2) does ...
            %
            % Inputs:
            %   arg1 - Description of first argument
            %   arg2 - Description of second argument
            %
            % Outputs:
            %   obj - The constructed ClassName object
        end
        
 	function output = publicMethod(obj, input)
            % PUBLICMETHOD Summary of public method purpose.
            %
            %   output = publicMethod(obj, input) does ...
            %
            % Inputs:
            %   obj - The ClassName object
            %   input - Description of input
            %
            % Outputs:
            %   output - Description of output
            %
            % Examples:
            %   obj = publicMethod(arg1, arg2);
            %   result = obj.method1(input);
        end
    end
    methods (Access = protected)
    	function output = protectedMethod(obj, input)
    	   % PROTECTEDMETHOD Summary of protected method purpose.
            %
            %   output = protectedMethod(obj, input) does ...
            %
            % Inputs:
            %   obj - The ClassName object
            %   input - Description of input
            %
            % Outputs:
            %   output - Description of output
    	end
    end
    methods (Access = private)
    	function output = privateMethod(obj, input)
    		% PRIVATEMETHOD Summary of private method purpose.
    	end
    end
end
```

