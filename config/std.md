# MATLAB Coding Standard  

Writing MATLAB code with consistent style, naming, and documentation.  

---

## 1. General Principles for Docstrings  

* **Public functions/methods**: Include full docstring with a complete description.  
* **Protected/private methods**: Only a **one-line description** is needed.  
* **Dependent property setters/getters**: Docstring should be **one line** and follow the format `SET.Property` / `GET.Property`.  
* Notes and See also sections are optional.  
* Capitalize the function/class name in the first line.  
* If you mention an input or output in the description, prefix it with @a to indicate it is an argument or return value.  
* If a comment or docstring line exceeds 80 characters, break it into multiple lines at logical points.  
* In-function comments are not allowed unless it is necessary.  
* **Math formulas in comments or docstrings**: Use LaTeX-like format; `\f$ ... \f$` for inline math, `\n\f[\n ... \n\f]\n` for display math.  

---

## 2. Function Docstring Template  

```matlab
function [output1, output2] = functionName(input1, input2, varargin)
% FUNCTIONNAME One-line description of function purpose.
%
%   [output1, output2] = functionName(input1, input2) performs a 
%   specific computation or operation.
%
%   [output1, output2] = functionName(input1, input2, option1) performs 
%   another computation or operation.
%
%   Any special considerations, limitations, or implementation details.
%
%   Detailed description should explain the behavior, the steps
%   performed, and how it uses @a input1 and @a input2 to produce @a output1
%   and @a output2. 
%
%   Example math: \f$ area = pi * r^2 \f$ or display math:
%
%   \f[
%       E = mc^2
%   \f]
%
% See also:
%   otherFunction1, otherFunction2

end
```

### Function Naming Conventions  

* Use **verbs** describing action: `computeArea`, `normalizeData`.  
* Use **camelCase**, must match the file name.  
* Avoid single-letter names except in short loops or formulas.  

---

## 3. Class Docstring Template  

```matlab
classdef ClassName < ParentClass
    % CLASSNAME One-line description of class purpose.
    %
    %   ClassName represents [short description of functionality].
    %   The detailed description should explain major capabilities,
    %   usage patterns, and key behaviors.
    %
    %   Any special considerations, limitations, or implementation details.
    %
    % See also:
    %   RelatedClass1, relatedFunction1

    properties
        Property   % Regular property, use PascalCase
    end
    
    properties (Dependent)
        DepProperty    % Dependent property, computed from others
    end

    methods
        function obj = ClassName(arg1, arg2)
            % CLASSNAME Construct an instance of ClassName with a complete description.
            %
            %   obj = ClassName(arg1) creates an object with @a arg1. 
        end

        function output = publicMethod(obj, input)
            % PUBLICMETHOD Performs a specific operation and returns the result.
            %
            %   output = publicMethod(input) performs a specific computation or operation.
        end
    end

    methods
        function obj = set.DepProperty(obj, value)
            % SET.DepProperty Sets the value of the dependent property 'DepProperty'.
        end

        function value = get.DepProperty(obj)
            % GET.DepProperty Returns the value of the dependent property 'DepProperty'.
        end
    end

    methods (Access = protected)
        function output = protectedMethod(obj, input)
            % PROTECTEDMETHOD One-line description only.
        end
    end

    methods (Access = private)
        function output = privateMethod(obj, input)
            % PRIVATEMETHOD One-line description only.
        end
    end
end
```

### Class Naming Conventions  

* Class names: **PascalCase** (`CircleManager`, `DataAnalyzer`).  
* Properties: **PascalCase** (`Radius`, `NDims`, `VariableNames`).  
* Methods: **camelCase verbs** (`getArea`, `plotCircle`).  
* Protected/private methods: one-line description only.  
* Getter docstrings: `GET.Property`, one-line description.  
* Setter docstrings: `SET.Property`, one-line description.  

---

## 4. Constants / Enumerations Naming  

* Prefer **PascalCase**, consistent with MATLAB toolboxes (`MaxEpochs`, `ExecutionEnvironment`).  
* UPPERCASE may be used for mathematical/physical constants (`PI`, `MAX_ITERATIONS`) if clarity is improved.  

---

## 5. Code Style  

* Indentation: 4 spaces per level, no tabs.  
* Line breaks: maximum 80 characters, break lines at logical points.  
* Operator spacing: spaces around operators (`=`, `+`, `-`, etc.).  
* Logical expressions: use consistent style (`&&` vs `&`, `||` vs `|`).  

---

## 6. Vectorization  

* Prefer vectorized operations over loops for performance.  
* Avoid growing arrays inside loops; preallocate memory when necessary.  
* Use built-in MATLAB functions for computations when possible.  

---

## 7. Argument Validation  

* For **all public methods and functions (except property getters)**,  
    use an `arguments` block for input validation and default values.  

* Do not use manual `nargin` checks or `inputParser` unless unavoidable.  

Example:  

```matlab
function result = computeNorm(x, p)
% COMPUTENORM Compute the p-norm of vector x.

    arguments
        x double {mustBeVector}
        p double {mustBePositive} = 2
    end

    result = sum(abs(x).^p)^(1/p);
end
```

---

## 8. Function Call Syntax  

* Always use **`name=value` syntax** for optional arguments when calling functions.  
* Avoid positional arguments for optional parameters.  
* This ensures readability, clarity, and forward compatibility.  

Example:  

```matlab
% Good (modern MATLAB style)
net = trainNetwork(XTrain, YTrain, layers, ...
    Options=trainingOptions("sgdm", ...
        MaxEpochs=20, ...
        InitialLearnRate=0.01));

% Bad (classical positional style)
net = trainNetwork(XTrain, YTrain, layers, trainingOptions("sgdm", 20, 0.01));
```

---

## 9. Summary Table of Conventions  

| Type                      | Convention                                                   |
| ------------------------- | ------------------------------------------------------------ |
| Function                  | Full docstring, **camelCase verbs**                          |
| Public method             | Full docstring, **camelCase verbs**, use `arguments` block   |
| Protected method          | One-line docstring, **camelCase method name**                |
| Private method            | One-line docstring, **camelCase method name**                |
| Dependent property getter | One-line docstring, format `GET.Property`, **PascalCase property** |
| Dependent property setter | One-line docstring, format `SET.Property`, **PascalCase property** |
| Regular properties        | **PascalCase** (e.g., `NDims`, `VariableNames`)              |
| Constants / Enumerations  | Prefer **PascalCase**, UPPERCASE optional for math/physics constants |
| Class                     | **PascalCase class name**                                    |
| Variable                  | **camelCase name**                                           |
| Optional args             | Always use **`name=value` syntax** in calls                  |
| Comment lines             | Max 80 characters, wrap sensibly                             |
| In-function comments      | Avoid unless absolutely necessary                            |
| Math formulas             | Use LaTeX-like format: `\f$ ... \f$` inline, `\n\[\n ... \n\]\n` for display math |

---

*This standard aligns with modern MATLAB toolbox practices and emphasizes clarity, consistency, and maintainability.*  
