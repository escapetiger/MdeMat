function values = referenceValues(space, reference, nodes, time)
%REFERENCEVALUES Evaluate a declared spatial reference at physical FE points.
%   A reference is numeric, static f(x), f(x,t), or f(x,varargin) with named
%   t. The returned layout is components-by-points. This performs no solve or
%   coefficient fitting and leaves kinetic reconstruction to its owner.
arguments
    space (1, 1) mdemat.common.fem.FiniteElementSpace
    reference
    nodes {mustBeNumeric}
    time = []
end
callback = reference;
args = {};
if isa(reference, 'function_handle') && ~isempty(time)
    arity = nargin(reference);
    if arity == 2
        callback = @(x) reference(x, time);
    elseif arity ~= 1
        args = {'t', time};
    end
end
values = space.feval(callback, nodes, args=args).';
end
