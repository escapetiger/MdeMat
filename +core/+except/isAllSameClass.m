function tf = isAllSameClass(varargin)
% ISALLSAMECLASS Check if all inputs belong to the same class.
%
%   tf = isAllSameClass(X1, X2, ...) returns true if all inputs belong to
%   the same class. For a single input argument, always returns true.
%
% See also:
%   core.except.isAllClass

firstClass = class(varargin{1});
tf = all(cellfun(@(x) strcmp(class(x), firstClass), varargin));

end