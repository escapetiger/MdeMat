function runToolsDesignReview(root, output)
%RUNTOOLSDESIGNREVIEW Observe current contracts without changing repository source.
oldPath = path; pathGuard = onCleanup(@() path(oldPath)); %#ok<NASGU>
addpath(root); setupMdeMat();
oldVisible = get(groot,'DefaultFigureVisible');
oldCurrent = get(groot,'CurrentFigure'); oldFigures = findall(groot,'Type','figure');
graphicsGuard = onCleanup(@() restoreGraphics(oldVisible,oldCurrent,oldFigures)); %#ok<NASGU>
set(groot,'DefaultFigureVisible','off');
work = tempname; mkdir(work); workGuard = onCleanup(@() rmdir(work,'s')); %#ok<NASGU>
names = {'config_partial_update','config_integrator_divergence','config_paths_collision', ...
    'analysis_selected_sequence','fem_unknown_reduction','graph_isolated_subsample', ...
    'visual_field_coordinate_collision','checkpoint_nested_exact','checkpoint_observer_identity', ...
    'rates_adjacent_large_scales','vlasov_ignored_exact','container_handle_alias', ...
    'external_method_inheritance','checkpoint_filename_collision','mixed_norm_definitions'};
operations = {@configPartial,@configStale,@() pathsCollision(work),@sequenceSelection, ...
    @unknownReduction,@isolatedGraph,@reservedField,@nestedExact,@observerIdentity, ...
    @closeRates,@vlasovExact,@handleAlias,@() externalMethod(work),@checkpointName,@mixedNorms};
probes = repmat(struct('Name','','Executed',false,'Observation',struct(), ...
    'HarnessError','','Diagnostic',''),1,numel(names));
for index = 1:numel(names)
    probes(index).Name = names{index};
    try
        probes(index).Observation = operations{index}();
        probes(index).Executed = true;
    catch failure
        probes(index).HarnessError = failure.identifier;
        probes(index).Diagnostic = getReport(failure,'extended','hyperlinks','off');
    end
end
report = struct('Scope','Design-review observations and contract demonstrations; no source repairs.', ...
    'MatlabRelease',version('-release'),'Computer',computer, ...
    'CapturedUTC',char(datetime('now',TimeZone='UTC')),'Probes',probes);
fid=fopen(output,'wt','n','UTF-8'); assert(fid>=0); closeGuard=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
fprintf('DESIGN_PROBES executed=%d harnessErrors=%d\n',sum([probes.Executed]),sum(~[probes.Executed]));
end
function result = configPartial()
s = mdemat.internal.execution.Scheme(); s.addConfig('a',default=1,validator=@(x) x>0);
s.addConfig('b',default=1,validator=@(x) x>0);
identifier='';try,s.loadFromStruct(struct('a',2,'b',-1));catch failure,identifier=failure.identifier;end
result=struct('RejectedIdentifier',identifier,'FirstValueAfterFailure',s.Config.a, ...
    'SecondValueAfterFailure',s.Config.b,'ExpectedUncommittedFirstValue',1);
end
function result = configStale()
s=mdemat.internal.execution.OdeScheme(config=struct('tFinal',1,'tOdeIntName','be'));
s.load(config=struct('tFinal',2,'tOdeIntName','fe'));
result=struct('ConfiguredFinal',s.Config.tFinal,'RuntimeFinal',s.TDisc.Timeline.Final, ...
    'ConfiguredMethod',{s.Config.tOdeIntName},'RuntimeClass',class(s.TDisc));
end
function result = pathsCollision(work)
file=fullfile(work,'paths.txt');writeText(file,sprintf('OPTION = paths = 7\nMATLABPATH = "%s"\n',work));
identifier='';try,value=mdemat.tools.environment.ConfigParser.parseFile(file);catch failure,identifier=failure.identifier;value=struct();end
result=struct('Identifier',identifier,'Value',value, ...
    'Note','OPTION paths and parser-produced resource paths occupy the same key namespace.');
end
function result = sequenceSelection()
a=mdemat.tools.analysis.Analyzer(nLevels=2);a.setComponents(struct('U',[2,4]));
a.setReductions(struct('U',struct('x',{{'L2'}},'v','sequence')));a.initialize();
a.Resolutions=[1;.5];a.Errors.U={sqrt(5);sqrt(5)/4};table=a.analyze();
identifier='';try,a.draw();catch failure,identifier=failure.identifier;end
result=struct('ReportedCount',a.NComponents,'ResultNames',{fieldnames(table)}, ...
    'DrawIdentifier',identifier,'ExpectedAggregateCount',1);
end
function result = unknownReduction()
space=makeSpace();a=mdemat.common.fem.FiniteElementAbsoluteMetric(space,reductions={'TYPO'});
a.setPoints(4);count=size(space.eval(a.Nodes,ones(space.NGlobalDofs,1)),1);
U=ones(count,1);V=zeros(count,1);identifier='';
try,value=a.eval(U,V);catch failure,identifier=failure.identifier;value=[];end
result=struct('Accepted',isempty(identifier),'Identifier',identifier, ...
    'ReturnedShape',size(value),'ReturnedNumel',numel(value));
end
function space = makeSpace()
space=mdemat.common.fem.FiniteElementSpace(mdemat.data.mesh.UniformGrid([0,1],2), ...
    mdemat.common.fem.BH1OrthotopeElement.modal(1,2));
end
function result = isolatedGraph()
g=mdemat.data.mesh.SeparableGraph({1:200},{zeros(0,2)});
p=mdemat.tools.visual.SeparableGraphPlotter();identifier='';
try,p.plot(g,MaxVertices=12,ShowEdges=false);catch failure,identifier=failure.identifier;end
result=struct('Identifier',identifier,'Vertices',g.NTotalVertices,'Edges',g.NTotalEdges);
end
function result = reservedField()
v=mdemat.tools.visual.Visualizer();v.setDensity(1);v.setComponents(struct('x',1));
v.addDataset('A','numeric');v.setRefXNodes(mdemat.data.geometry.Orthotope([0,1]));
v.setPhyXNodes(mdemat.data.mesh.UniformGrid([0,1],2));ds=v.Database.getDataset('A');before=ds.Data.x;
sampler=struct('eval',@(~,u)u,'feval',@(varargin)[]);v.update(sampler,struct('x',[5;6]));
result=struct('CoordinatesBefore',before,'SameSlotAfterFieldUpdate',ds.Data.x, ...
    'CoordinateSlotWasOverwritten',~isequal(before,ds.Data.x));v.finalize();
end
function result = nestedExact()
identifier='';try,mdemat.internal.execution.checkpointConfiguration(struct('exact',struct('F',@sin)));catch failure,identifier=failure.identifier;end
result=struct('Identifier',identifier,'Note','This nested exact form is admitted by RadiationScheme.');
end
function result = observerIdentity()
a=mdemat.internal.execution.checkpointConfiguration(struct('dt',.1,'components',struct('U',1)));
b=mdemat.internal.execution.checkpointConfiguration(struct('dt',.1,'components',struct('U',2)));
result=struct('SameNumericalDt',true,'IdentityEqual',isequal(a,b), ...
    'ObservationFieldsIncluded',{fieldnames(a.values)});
end
function result = closeRates()
h=[1e200;1e200*(1-eps)];a=mdemat.tools.analysis.convergenceRates(h,[2;1]);
expected=log(2)/log1p((h(1)-h(2))/h(2));
result=struct('StrictlyDecreasing',all(diff(h)<0),'LogDifference',log(h(1))-log(h(2)), ...
    'ActualRate',sprintf('%.17g',a.LocalRates),'Defined',a.Defined, ...
    'FiniteLog1pReference',expected);
end
function result = vlasovExact()
a=mdemat.tools.analysis.Analyzer();
s=mdemat.application.plasma.VlasovScheme(config=struct('exact',@sin,'analyzer',a));
result=struct('ConfiguredExactIsFunction',isa(s.Config.exact,'function_handle'), ...
    'RegisteredReferenceFields',{fieldnames(a.Exacts)},'HasExact',a.HasExact);
end
function result = handleAlias()
a=mdemat.tools.visual.Dataset(data=struct('U',1));collection=mdemat.tools.visual.Database();
collection.setDataset('A',a);returned=collection.getDataset('A');returned.Data.U=9;
t=mdemat.tools.timing.Timer(verbose=false);value=struct('Timer',t);copy=value;copy.Timer.Verbose=true;
result=struct('DatabaseChangedThroughReturnedHandle',collection.getDataset('A').Data.U, ...
    'StructCopyStillSharesEmbeddedHandle',value.Timer.Verbose, ...
    'Classification','Ownership demonstration; a value struct is not a recursive deep copy.');
end
function result = externalMethod(work)
namespace=fullfile(work,'+audittoolsdesign');base=fullfile(namespace,'@PlotBase');mkdir(base);
writeText(fullfile(base,'PlotBase.m'),sprintf('classdef PlotBase\nmethods\nh=plot(obj,ax)\nend\nend\n'));
writeText(fullfile(base,'plot.m'),sprintf('function h=plot(~,ax)\nh=line(ax,[0,1],[0,1]);\nend\n'));
writeText(fullfile(namespace,'PlotChild.m'),sprintf('classdef PlotChild < audittoolsdesign.PlotBase\nend\n'));
old=path;guard=onCleanup(@()path(old)); %#ok<NASGU>
addpath(work);obj=audittoolsdesign.PlotChild();f=figure('Visible','off');ax=axes(f);h=obj.plot(ax);
result=struct('ExternalMethodRecognized',isgraphics(h,'line'), ...
    'InheritedBySubclass',isa(obj,'audittoolsdesign.PlotBase'),'ExplicitAxesRespected',h.Parent==ax);
delete(f);
end
function result = checkpointName()
times=[.11,.14];names=arrayfun(@(t)strrep(sprintf('t%.1f',t),'.','p'),times,UniformOutput=false);
result=struct('AcceptedTimes',times,'Names',{names},'NamesUnique',numel(unique(names))==numel(times), ...
    'Scope','Reproduces the filename conversion only; no Gaussian simulation or real checkpoint was written.');
end
function result = mixedNorms()
U=[1,0;0,1];w=[.5,.5];r=mdemat.tools.analysis.weightedError(U,zeros(2),w);
a=sum(abs(U).*w,2);b=max(abs(U),[],2);
result=struct('WeightedGlobalL1',r.L1,'LegacySequenceL1',norm(a,2), ...
    'WeightedGlobalLInf',r.LInf,'LegacySequenceMax',norm(b,2), ...
    'Classification','Distinct norm definitions; this does not by itself establish an existing calculation bug.');
end
function writeText(file,value)
f=fopen(file,'wt','n','UTF-8');assert(f>=0);guard=onCleanup(@()fclose(f)); %#ok<NASGU>
fprintf(f,'%s',value);
end
function restoreGraphics(visible,current,previous)
all=findall(groot,'Type','figure');for i=1:numel(all),if ~any(all(i)==previous),delete(all(i));end,end
set(groot,'DefaultFigureVisible',visible);if ~isempty(current)&&isgraphics(current),set(groot,'CurrentFigure',current);end
end
