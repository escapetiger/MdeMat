function runToolsReview(root, output)
oldPath = path;
pathGuard = onCleanup(@() path(oldPath)); %#ok<NASGU>
addpath(root); setupMdeMat();
visibility = get(groot, 'DefaultFigureVisible');
figures = findall(groot, 'Type', 'figure');
figureGuard = onCleanup(@() restoreFigures(figures, visibility)); %#ok<NASGU>
set(groot, 'DefaultFigureVisible', 'off');
work = tempname; mkdir(work);
folderGuard = onCleanup(@() rmdir(work, 's')); %#ok<NASGU>
names = {'predicate_multidimensional_cell','field_string_array', ...
    'field_matrix_output_arity','field_empty_name','same_class_empty', ...
    'assert_mixed_array','timer_matrix_names','timer_anonymous_collision', ...
    'hash_standard_and_collision','analysis_rate_overflow', ...
    'weighted_error_scaled_extremes','weighted_relative_overflow_policy', ...
    'database_coordinate_order','visual_component_identity', ...
    'visual_history_nonuniform','visual_finalize_reuse','visual_render_failure', ...
    'figure_close_veto','config_power_association','config_source_relative_path', ...
    'io_complex_corrupt_bytecount'};
operations = {@predicate,@fieldsStrings,@fieldsMatrix,@fieldsEmpty, ...
    @sameClass,@assertMixed,@timerMatrix,@timerCollision,@hashProbe,@rateProbe, ...
    @weightedProbe,@relativeProbe,@coordinateProbe,@componentProbe, ...
    @historyProbe,@reuseProbe,@renderProbe,@closeProbe, ...
    @() powerProbe(work),@() configPathProbe(work),@() ioProbe(work)};
records = repmat(struct('name','','executed',false,'observation',struct(), ...
    'errorIdentifier','','diagnostic',''),1,numel(names));
for index = 1:numel(names)
    records(index).name = names{index};
    try
        records(index).observation = operations{index}();
        records(index).executed = true;
    catch exception
        records(index).errorIdentifier = exception.identifier;
        records(index).diagnostic = getReport(exception,'extended','hyperlinks','off');
    end
end
report = struct('scope','Read-only planning audit probes; findings are observations, not source repairs.', ...
    'release',version('-release'),'computer',computer, ...
    'capturedUTC',char(datetime('now',TimeZone='UTC')),'probes',records);
writeJson(output,report);
fprintf('AUDIT_PROBES executed=%d unexpectedProbeErrors=%d\n',sum([records.executed]),sum(~[records.executed]));
end

function r = predicate()
r = struct('actual',mdemat.tools.base.isAllClass({1,2;3,'bad'},'double'),'expected',false);
end
function r = fieldsStrings()
r = struct('actual',mdemat.tools.base.hasFields(struct('a',1,'b',2),["a","b"]),'expected',true);
end
function r = fieldsMatrix()
s = struct('a',1,'b',2,'c',3); f = {'a','b';'c','missing'};
one = mdemat.tools.base.hasFields(s,f);
[two,missing] = mdemat.tools.base.hasFields(s,f);
r = struct('singleOutput',one,'doubleOutput',two,'missing',missing,'expectedScalar',false);
end
function r = fieldsEmpty()
one = mdemat.tools.base.hasFields(struct('a',1),{''});
[two,missing] = mdemat.tools.base.hasFields(struct('a',1),{''});
r = struct('singleOutput',one,'doubleOutput',two,'missing',missing,'expected',false);
end
function r = sameClass()
try
    actual = mdemat.tools.base.isAllSameClass();
    r = struct('accepted',true,'actual',actual,'errorIdentifier','');
catch e
    r = struct('accepted',false,'actual',[],'errorIdentifier',e.identifier);
end
end
function r = assertMixed()
try
    mdemat.tools.except.assert([true,false],'AuditMixed','Must not pass mixed values.');
    r = struct('threw',false,'identifier','');
catch e
    r = struct('threw',true,'identifier',e.identifier);
end
end
function r = timerMatrix()
t = mdemat.tools.timing.Timer(names={'a','b';'c','d'},verbose=false);
r = struct('actualCount',numel(t.Records),'expectedCount',4);
end
function r = timerCollision()
t = mdemat.tools.timing.Timer(verbose=false);
t.start(record='Record_2'); t.stop(); t.start();
names = {t.Records.Name};
try
    t.stop(); r = struct('names',{names},'stopSucceeded',true,'errorIdentifier','');
catch e
    r = struct('names',{names},'stopSucceeded',false,'errorIdentifier',e.identifier);
end
end
function r = hashProbe()
a = uint8([127,128,255]); b = uint8([127,127,127]);
raw = mdemat.tools.base.hash(a,format='uint8');
other = mdemat.tools.base.hash(b,format='uint8');
canonical = lower(reshape(dec2hex(raw,2).',1,[]));
r = struct('header','uint8|1,3,|sparse=0|complex=0|', ...
    'payload',double(a),'hex',mdemat.tools.base.hash(a), ...
    'rawHex',canonical,'collisionWithSaturatedPayload',isequal(raw,other));
end
function r = rateProbe()
a = mdemat.tools.analysis.Analyzer(nLevels=2);
a.setComponents(struct('U',1)); a.initialize(); a.Resolutions = [1;0.1];
a.Errors.U = {1e200*ones(3,1);1e-200*ones(3,1)};
t = a.analyze();
stable = mdemat.tools.analysis.convergenceRates([1;0.1],[1e200;1e-200]);
r = struct('analyzerRate',char(string(t.U1.L2_Order(2))), ...
    'standaloneRate',stable.LocalRates,'expected',400);
end
function r = weightedProbe()
w = 1e-308;
a = mdemat.tools.analysis.weightedError(realmax,-realmax,w);
expected = 2*(realmax*sqrt(w));
r = struct('actualL2',sprintf('%.17g',a.L2),'finiteWeightedL2Reference',expected, ...
    'actualRelative',sprintf('%.17g',a.RelativeL2),'expectedRelative',2);
end
function r = relativeProbe()
a = mdemat.tools.analysis.weightedError(1,realmin/16,1);
r = struct('relative',sprintf('%.17g',a.RelativeL2),'relativeDefined',a.RelativeDefined, ...
    'finiteNorms',a.FiniteNorms,'classification','Contract decision: ratio representability versus mathematical definition.');
end
function r = coordinateProbe()
d = mdemat.tools.visual.Database(); x = [0,1,0,1;1,1,0,0];
ds = mdemat.tools.visual.Dataset(type='fixed',data=struct('x',x,'U',(x(1,:)+10*x(2,:)).'));
d.setDataset('A',ds); grouped = d.groupBy('U1'); result = grouped.getDataset('A');
r = struct('actual',result.Data.u,'expected',[0,10;1,11],'axes',{result.Data.x});
end
function r = componentProbe()
v = mdemat.tools.visual.Visualizer();
v.setDensity(1); v.setComponents(struct('u',[2,4]));
v.addDataset('A','numeric',struct('x',0));
v.setRefXNodes(mdemat.data.geometry.Orthotope([-1,1]));
sampler = struct('eval',@(~,u) u,'feval',@(varargin) []);
v.update(sampler,struct('u',[10,20,30,40]));
two = v.Database.groupBy('u2'); ds = two.getDataset('A'); four = v.Database.groupBy('u4');
r = struct('selectedTwo',ds.Data.u,'expectedTwo',20, ...
    'selectedFourDatasetCount',four.NDatasets,'expectedCount',1);
v.finalize();
end
function r = historyProbe()
v = mdemat.tools.visual.Visualizer(); v.setDensity(1);v.setComponents(struct('u',1));
v.addDataset('A','numeric');v.setRefXNodes(mdemat.data.geometry.Orthotope([-1,1]));
v.setTimeline(.3,3);v.Timeline.Now=.2;
t=mdemat.data.mesh.DynamicTimeline(.3,nSteps=2);t.StepSizeQueue=[.2,.1,0];t.Next=.3;
sampler=struct('eval',@(~,u) u,'feval',@(varargin) []);
v.update(sampler,struct('u',.3),timeline=t,V0={struct('u',.1),struct('u',0)},n0=2);
ds=v.Database.getDataset('A'); r=struct('actual',ds.Data.u,'expected',.2);
v.finalize();
end
function r = reuseProbe()
v=mdemat.tools.visual.Visualizer();v.setDensity(1);v.setComponents(struct('U',1));
v.initialize();v.finalize();v.addDataset('A','numeric');v.addPlotter('p','1d');v.initialize();
r=struct('axesEmptyAfterReinitialize',isempty(v.Plotters.p.Axes),'status',double(v.Status));
v.finalize();
end
function r = renderProbe()
v=mdemat.tools.visual.Visualizer();v.setComponents(struct('U',1));v.addPlotter('p','1d');
id=''; try,v.render();catch e,id=e.identifier;end
r=struct('errorIdentifier',id,'renderingFlagAfterFailure',bitand(v.Status,v.IsRendering)~=0);
v.finalize();
end
function r = closeProbe()
m=mdemat.tools.visual.FigureManager();h=m.createFigures(1);f=h{1};
guard=onCleanup(@() delete(f)); %#ok<NASGU>
set(f,'CloseRequestFcn',@(~,~) []);m.closeAll();
r=struct('stillAlive',isgraphics(f),'managerCount',m.NFigures);
end
function r = powerProbe(work)
file=fullfile(work,'power.txt');writeText(file,'OPTION = power = 2^3^2');
parsed=mdemat.tools.environment.ConfigParser.parseFile(file);
r=struct('parser',parsed.power,'matlab',2^3^2);
end
function r = configPathProbe(work)
a=fullfile(work,'source');b=fullfile(work,'caller');mkdir(a);mkdir(b);
mkdir(fullfile(a,'shared'));mkdir(fullfile(b,'shared'));
file=fullfile(a,'config.txt');writeText(file,'MATLABPATH = shared');
old=pwd;restore=onCleanup(@() cd(old)); %#ok<NASGU>
cd(b);one=mdemat.tools.environment.ConfigParser.parseFile(file);
value=one.paths{1};if strcmp(value,'shared'),first=fullfile(b,value);else,first=value;end
cd(work);two=mdemat.tools.environment.ConfigParser.parseFile(file);
r=struct('firstBinding',first,'secondBinding',two.paths{1},'expected',fullfile(a,'shared'));
end
function r = ioProbe(work)
file=fullfile(work,'complex.json');
mdemat.tools.io.writeRecord(file,[1+2i,3+4i]);
record=jsondecode(fileread(file));names=fieldnames(record.arrays);name=names{1};
record.arrays.(name).imaginary=matlab.net.base64encode(typecast(double(9),'uint8'));
writeJson(file,record);
try
    values=mdemat.tools.io.readRecord(file);
    r=struct('corruptRecordAccepted',true,'real',real(values),'imaginary',imag(values),'errorIdentifier','');
catch e
    r=struct('corruptRecordAccepted',false,'real',[],'imaginary',[],'errorIdentifier',e.identifier);
end
end
function writeText(file,text)
h=fopen(file,'w','n','UTF-8');guard=onCleanup(@() fclose(h)); %#ok<NASGU>
fprintf(h,'%s\n',text);
end
function writeJson(file,data)
writeText(file,jsonencode(data,PrettyPrint=true));
end
function restoreFigures(original,visibility)
delete(setdiff(findall(groot,'Type','figure'),original));
set(groot,'DefaultFigureVisible',visibility);
end
