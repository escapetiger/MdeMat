function tf = isAllSameClass(varargin)
%ISALLSAMECLASS Check exact class equality, without subclass matching.
%
%   tf = isAllSameClass(X1, X2, ...) returns true if all inputs belong to
%   the same class. Zero or one input returns true. An empty array still has
%   its declared class; this checks input classes, not array contents.
%
% See also:
%   mdemat.tools.base.isAllClass

if nargin == 0, tf = true; return; end
firstClass = class(varargin{1});
tf = all(cellfun(@(x) strcmp(class(x), firstClass), varargin));

end
