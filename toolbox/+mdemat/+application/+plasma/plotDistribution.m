function fields = plotDistribution(state, plan, components, velocityBounds, options)
%PLOTDISTRIBUTION Reconstruct physical f(x,v) as a scalar phase-space view.
%   F includes the model's Hermite basis and (for AP) sqrtRhoInf scaling.
%   Distance explicitly denotes f-f_equilibrium. Tools receive no State.
%   A supplied reference uses aligned Cartesian point pairs f(x,v), f(x,v,t)
%   or f(x,v,varargin) with named t, returning one scalar per pair.
arguments
    state mdemat.common.fem.KineticState
    plan (1,1) struct
    components (1,1) struct
    velocityBounds (1,2) double {mustBeFinite}
    options.reference = []
    options.time = []
end
fields=struct();
requestedF=isfield(components,'F')&&~isempty(components.F);
requestedDistance=isfield(components,'Distance')&&~isempty(components.Distance);
if ~requestedF && ~requestedDistance,return;end
if state.XDisc.NDims~=1 || ~isfield(plan,'Axes')
    error('mdemat:plasma:plotDistribution:Geometry','Phase-space display requires a 1D tensor spatial grid.');
end
velocity=linspace(velocityBounds(1),velocityBounds(2),512);
if isempty(options.reference)
    values=state.distribution(plan.Nodes,velocity);
    values=values(plan.Permutation,:);
else
    x=plan.Axes{1};[X,V]=ndgrid(x,velocity);
    raw=mdemat.application.plasma.referenceDistribution( ...
        options.reference,X(:).',V(:).',options.time);
    values=reshape(raw,size(X));
end
axes={plan.Axes{1},velocity};
if requestedF
    fields.F=struct('Values',reshape(values,1,[]),'ComponentIds',1,'Axes',{axes});
end
if requestedDistance
    if ~ismethod(state,'equilibrium')
        error('mdemat:plasma:plotDistribution:Equilibrium','This model does not supply an equilibrium.');
    end
    equilibrium=state.equilibrium(plan.Nodes,velocity);
    equilibrium=equilibrium(plan.Permutation,:);
    fields.Distance=struct('Values',reshape(values-equilibrium,1,[]), ...
        'ComponentIds',1,'Axes',{axes});
end
end
