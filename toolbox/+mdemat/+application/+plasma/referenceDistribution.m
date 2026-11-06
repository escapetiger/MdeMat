function values = referenceDistribution(reference, x, velocity, time)
%REFERENCEDISTRIBUTION Evaluate a physical scalar on aligned (x,v) point pairs.
%   F is numeric or f(x,v), f(x,v,t), or f(x,v,varargin) with named t.
%   A scalar output broadcasts; otherwise exactly one value per pair is needed.
%   This does not form a Cartesian grid, fit coefficients, or run a solver.
%   Output finiteness is checked by assessment; plotting may retain missing NaN.
arguments
    reference
    x (1,:) double {mustBeReal,mustBeFinite}
    velocity (1,:) double {mustBeReal,mustBeFinite}
    time {mustBeNumeric,mustBeScalarOrEmpty,mustBeReal,mustBeFinite} = []
end
if numel(x) ~= numel(velocity)
    error('mdemat:plasma:referenceDistribution:Shape','Physical point pairs must align.');
end
if isnumeric(reference)
    raw = reference;
elseif isa(reference,'function_handle')
    arity = nargin(reference);
    if arity == 2
        raw = reference(x,velocity);
    elseif arity == 3
        if isempty(time)
            error('mdemat:plasma:referenceDistribution:Time','The physical reference requires a time.');
        end
        raw = reference(x,velocity,time);
    elseif arity < 0
        raw = reference(x,velocity,'t',time);
    else
        error('mdemat:plasma:referenceDistribution:Arity','Use f(x,v), f(x,v,t), or named t.');
    end
else
    error('mdemat:plasma:referenceDistribution:Type','Reference must be numeric or a callback.');
end
if ~isnumeric(raw) || ~(isscalar(raw) || numel(raw) == numel(x))
    error('mdemat:plasma:referenceDistribution:Shape', ...
        'The physical reference must return one scalar per point pair.');
end
if isscalar(raw), raw = repmat(raw,1,numel(x)); end
values = double(raw(:).');
end
