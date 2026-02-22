function tf = isAllSameClass(varargin)
% ISALLSAMECLASS Check if all inputs belong to the same class.
%
%   tf = isAllSameClass(X1, X2, ...) returns true if all inputs belong to
%   the same class. For a single input argument, always returns true.
%
% Inputs:
%   varargin - Input arguments lists
%
% Outputs:
%   tf - True if all input arguments have the same class
%
% Examples:
%   % Check if all inputs are the same numeric type
%   result = isAllSameClass(1.0, 2.0, 3.0);  % returns true
%
%   % Check mixed types
%   result = isAllSameClass(1.0, int32(2), 3.0);  % returns false
%
%   % Check objects of the same class
%   obj1 = containers.Map();
%   obj2 = containers.Map();
%   result = isAllSameClass(obj1, obj2);  % returns true
%
%   % Single input always returns true
%   result = isAllSameClass(42);  % returns true
%
% See Also:
%   core.validate.isAllClass

firstClass = class(varargin{1});
tf = all(cellfun(@(x) strcmp(class(x), firstClass), varargin));

end