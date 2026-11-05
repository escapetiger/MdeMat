# MdeMat Development Guide

## 1. Project Overview

### 1.1 Purpose

MdeMat is a modular and scalable MATLAB library for modern researches on numerical algorithms for differential equations.

### 1.2 Design Philosophy

The project is founded on four core principles:

- **Modularity**: Decomposing complex computational tasks into manageable, reusable components
- **Layered Architecture**: Clear separation of concerns through dependency-based layers
- **Object-Oriented Design**: Leveraging MATLAB's object-oriented programming capabilities
- **Code Quality**: Adhering to rigorous software development standards

## 2. Project Architecture

MdeMat is organized into the following computational layers with strict upward-only dependencies:

```
Level 5 (Physics):        +physics/
Level 4 (Profilers):      +profilers/  
Level 3 (Discretization): +fem/
Level 2 (Approximation):  +approx/
Level 1 (Mathematics):    +core/ (Level 1 components)
Level 0 (Infrastructure): +core/ (Level 0 components)
```

The architecture enforces strict dependency rules to maintain modularity:

- **Level 5** (Physics: `+physics/`) can use: Levels 0, 1, 2, 3, 4
- **Level 4** (Profilers: `+profilers/`) can use: Levels 0, 1, 2, 3
- **Level 3** (Discretization: `+fem/`) can use: Levels 0, 1, 2
- **Level 2** (Approximation: `+approx/`) can use: Levels 0, 1
- **Level 1** (Core Math: `+core/` Level 1 components) can use: Level 0
- **Level 0** (Core Infrastructure: `+core/` Level 0 components): Minimal external dependencies

## 3. Coding Standards

### 3.1 Code Style Guidelines

All developers must adhere to the principles outlined in Richard K. Johnson's "The Elements of MATLAB Style":

#### Naming Conventions

- Use meaningful and descriptive names
- Employ lowerCamelCase for variables and functions
- Use UpperCamelCase for class names
- Avoid cryptic abbreviations

#### Function Design

- Write focused, single-responsibility functions
- Limit function complexity
- Provide clear input and output specifications
- Include comprehensive documentation

#### Documentation Requirements

- Write detailed header comments for all functions and classes
- Describe purpose, inputs, outputs, and usage
- Document any side effects or special considerations
- Keep documentation concise but complete

### 3.2 Object-Oriented Programming Principles

- Prioritize composition over inheritance
- Maintain small, focused classes
- Encapsulate data and behavior
- Use private properties when appropriate
- Design clear and minimal public interfaces

### 3.4 Documentation Template

#### Function Documentation Template

```matlab
function [output1, output2] = functionName(input1, input2, varargin)
% FUNCTIONNAME One-line summary of function purpose.
%
% Detailed description of what the function does, algorithmic approach,
% and important considerations.
%
% Syntax:
%   output = functionName(input)
%   [out1, out2] = functionName(in1, in2, ...)
%   output = functionName(___,'ParameterName',ParameterValue)
%
% Inputs:
%   in1 - Description of first input argument
%   in2 - Description of second input argument
%   'ParameterName' - Description of parameter
%
% Outputs:
%   out1 - Description of first output
%   out2 - Description of second output
%
% Examples:
%   % Basic usage
%   result = functionName(input);
%
%   % Advanced usage with parameters
%   [out1, out2] = functionName(in1, in2, 'Parameter', value);
%
% Notes:
%   Any special considerations, limitations, or implementation details.
%
% See Also:
% 	packageName.moduleName.relatedFunction1, packageName.moduleName.relatedFunction2
end
```

#### Class Documentation Template

```matlab
classdef ClassName < ParentClass
    % CLASSNAME One-line summary of class purpose.
    %
    % Extended description of class purpose and functionality.
    % This section should outline the major capabilities and use cases.
    %
    % Examples:
    %   obj = ClassName(arg1, arg2);
    %   result = obj.method1(input);
    %
    % See also:
    %   packageName.moduleName.relatedClass, packageName.moduleName.relatedFunction
    
    properties
        % Description of property
        PropertyName
        
        % Description with expected values or constraints
        AnotherProperty
    end
    
    methods
        function obj = ClassName(arg1, arg2)
            % CLASSNAME Constructor for ClassName
            %
            % Syntax:
            %   obj = ClassName(arg1, arg2)
            %
            % Inputs:
            %   arg1 - Description of first argument
            %   arg2 - Description of second argument
            %
            % Outputs:
            %   obj - The constructed ClassName object
        end
        
        function output = methodName(obj, input)
            % METHODNAME Summary of method purpose
            %
            % Syntax:
            %   output = methodName(obj, input)
            %
            % Inputs:
            %   obj - The ClassName object
            %   input - Description of input
            %
            % Outputs:
            %   output - Description of output
        end
    end
end
```

