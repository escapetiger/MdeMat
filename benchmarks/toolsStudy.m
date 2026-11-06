function study = toolsStudy(options)
%TOOLSSTUDY Measure validated tool operations inside this MATLAB process.
%   One untimed validation run precedes warmup+repetitions measured by benchmark.
%   Oracles run outside timing. Inputs use a private stream; global RNG/path/cwd
%   are unchanged. Temporary I/O and the owned graphics figure are cleaned up.
%   Native memory is sampled before/after each case, never called a peak or an
%   allocation count. No previous-code baseline or hardware threshold is used.
arguments
    options.warmup (1,1) double {mustBeInteger,mustBeNonnegative,mustBeFinite} = 1
    options.repetitions (1,1) double {mustBeInteger,mustBePositive,mustBeFinite} = 3
    options.denseSize (1,2) double {mustBeInteger,mustBePositive,mustBeFinite} = [256,256]
    options.spatialCells (1,1) double {mustBeInteger,mustBePositive,mustBeFinite} = 32
    options.graphEdges (1,1) double {mustBeInteger,mustBePositive,mustBeFinite} = 128
    options.seed (1,1) double {mustBeInteger,mustBeNonnegative,mustBeFinite} = 7
end
stream = RandStream('mt19937ar','Seed',options.seed);
dense = rand(stream,options.denseSize);
imaginary = rand(stream,options.denseSize);
n = max(options.denseSize);
sparseData = spdiags(rand(stream,n,3),-1:1,n,n);
payloads = {dense,sparseData,complex(dense,imaginary)};
labels = {'dense','sparse','complex'};
folder = tempname;
[ok,message] = mkdir(folder);
if ~ok, error('mdemat:benchmark:toolsStudy:Directory','%s',message); end
filesGuard = onCleanup(@() removeFolder(folder)); %#ok<NASGU>
identity = mdemat.tools.environment.versionInfo();
cases = {};

points = prod(options.denseSize);
view = struct('Values',rand(stream,4,points),'ComponentIds',[2,4,6,8], ...
    'Axes',{{linspace(0,1,points)}});
collection = struct('NUMERIC',struct('Type','numeric','Data',struct('u2',view)));
expected = view.Values(2,:);
record('collection-selection',@() mdemat.tools.visual.selectComponent(collection,'u2',4), ...
    @(result) isequal(result.NUMERIC.Data.u2.Values,expected), ...
    'Select an original component ID from an existing value collection');
collection = constructCollection(view);
original = collection;
record('collection-construction',@() constructCollection(view), ...
    @(result) isequaln(result,original), ...
    'Construct eight datasets by three fields using ordinary value views');
record('collection-copy-update',@() updateCollection(collection), ...
    @(result) result.D01.Data.u2.Values(2,1)==original.D01.Data.u2.Values(2,1)+1 && ...
    isequaln(collection,original) && isequaln(result.D02,original.D02), ...
    'Explicit struct copy and one numeric entry update; original collection retained');
space = mdemat.common.fem.FiniteElementSpace( ...
    mdemat.data.mesh.UniformGrid([-1,1],options.spatialCells), ...
    mdemat.common.fem.BH1OrthotopeElement.modal(1,3));
plan = mdemat.common.fem.samplingPlan(space,4);
polynomial = @(x) [1+x;2-x+3*x.^2];
coefficients = space.fit(polynomial);
reference = polynomial(plan.Coordinates);
record('fem-plan',@() mdemat.common.fem.samplingPlan(space,4), ...
    @(result) abs(sum(result.Weights)-2)<1e-12 && all(result.Weights>0), ...
    'Construct numerical quadrature plan, physical coordinates and Jacobians');
record('fem-samples',@() mdemat.common.fem.sampleField(space,coefficients,plan), ...
    @(result) isequal(size(result),size(reference)) && max(abs(result-reference),[],'all')<1e-12, ...
    'Evaluate a representable quadratic field at an existing plan');
record('weighted-error',@() mdemat.tools.analyze.weightedError( ...
    reference+[1;-2],reference,plan.Weights), ...
    @(result) abs(result.L1-6)<1e-12 && abs(result.L2-sqrt(10))<1e-12 && ...
    abs(result.LInf-2)<1e-12, 'Physical norms only; sampling and oracle excluded');

for kind = 1:numel(payloads)
    payload = payloads{kind};
    for extension = {'mat','json'}
        file = fullfile(folder,[labels{kind},'.',extension{1}]);
        record([labels{kind},'-',extension{1}],@() roundtrip(file,payload), ...
            @(result) isequaln(result.Value,payload) && ...
            strcmp(class(result.Value),class(payload)) && ...
            issparse(result.Value)==issparse(payload) && isreal(result.Value)==isreal(payload), ...
            'Write/read complete record; explicit overwrite inside an owned temporary directory');
    end
end
smallLeaves = num2cell(1:64);
for extension = {'mat','json'}
    file = fullfile(folder,['small-leaves.',extension{1}]);
    record(['small-leaves-',extension{1}],@() roundtrip(file,smallLeaves), ...
        @(result) isequaln(result.Value,smallLeaves) && result.LeafCount==64, ...
        'Write/read sixty-four numeric scalar leaves; explicit overwrite in an owned directory');
end
A = spdiags(repmat([-1,3,-1],n,1),-1:1,n,n);
b = A*ones(n,1);
for enabled = [false,true]
    record(sprintf('linear-solve-recordHash-%d',enabled),@() solveOnce(A,b,enabled), ...
        @(result) norm(result.Value-1,Inf)<1e-11 && result.HashRecorded==enabled, ...
        'Construct, factor and solve the same SPD system; optional diagnostic hash');
end
record('timer-intervals',@timerIntervals,@validTimer, ...
    'Construct Timer and commit thirty named no-work intervals');
record('metadata-version',@mdemat.tools.environment.versionInfo, ...
    @(result) isequal(result,identity),'Read and validate installed identity metadata');
record('metadata-environment',@mdemat.tools.environment.describe, ...
    @(result) strcmp(result.ToolboxVersion,identity.ToolboxVersion) && ...
    ~result.InstalledProductsIncluded,'Describe native environment without product inventory');

previousFigure = get(groot,'CurrentFigure');
ownedFigure = figure('Visible','off');
graphicsGuard = onCleanup(@() removeFigure(ownedFigure,previousFigure)); %#ok<NASGU>
axis = axes('Parent',ownedFigure,'NextPlot','add','XLim',[-1.2,1.2],'YLim',[-1.2,1.2]);
angles = linspace(0,2*pi,options.graphEdges+1).';
vertices = [cos(angles),sin(angles)];
edges = [(1:options.graphEdges).',(2:options.graphEdges+1).'];
expectedEdges = [vertices(edges(:,1),:),vertices(edges(:,2),:)];
for batch = [false,true]
    record(sprintf('graph-batch-%d',batch),@() drawGraph(axis,vertices,edges,batch), ...
        @(result) isequal(result.Endpoints,expectedEdges) && ...
        result.ArtistCount==(double(batch)+(~batch)*size(edges,1)), ...
        'Create edge artists, drawnow, capture endpoints/count and delete owned artists; axes reused');
end
study = struct('Schema','mdemat.toolsStudy','SchemaVersion',1,'Options',options, ...
    'Inputs',struct('Dense',description(dense),'Sparse',description(sparseData), ...
    'Complex',description(payloads{3}),'SpatialQuadraturePoints',numel(plan.Weights), ...
    'GraphVertices',size(vertices,1),'GraphEdges',size(edges,1),'SmallLeafCount',64, ...
    'CollectionDatasets',8,'CollectionFields',3), ...
    'Cases',[cases{:}],'Environment',mdemat.tools.environment.describe(), ...
    'ExecutionStatus','COMPLETED','Validation','Every case oracle checked before and after timing', ...
    'MeasurementScope','Current MATLAB process; one untimed oracle run plus requested warmup/repetitions', ...
    'MemoryScope','Before/after process snapshots and their difference; not peak or allocation bytes', ...
    'PerformanceClaims','No old-code baseline or hardware acceptance threshold; no superiority inferred', ...
    'Deferred',{{'High-dimensional/low-rank slice memory: no new materialization/cache optimization; measure a representative application before changing it.'}});
    function record(name,operation,verify,scope)
        result = operation();
        assert(verify(result),'mdemat:benchmark:toolsStudy:Oracle','Case %s failed its oracle.',name);
        before = mdemat.tools.analyze.memorySnapshot();
        timing = mdemat.tools.analyze.benchmark(@execute,warmup=options.warmup, ...
            repetitions=options.repetitions,scope=scope);
        after = mdemat.tools.analyze.memorySnapshot();
        assert(verify(result),'mdemat:benchmark:toolsStudy:Oracle','Case %s changed its answer.',name);
        bytes = NaN; artists = NaN; leaves = NaN;
        if isstruct(result) && isfield(result,'SerializedBytes'), bytes = result.SerializedBytes; end
        if isstruct(result) && isfield(result,'ArtistCount'), artists = result.ArtistCount; end
        if isstruct(result) && isfield(result,'LeafCount'), leaves = result.LeafCount; end
        cases{end+1} = struct('Name',name,'Timing',timing,'Validated',true, ...
            'UntimedExecutions',1+options.warmup,'MemoryBefore',before,'MemoryAfter',after, ...
            'SnapshotDeltaBytes',after.Bytes-before.Bytes,'SerializedBytes',bytes, ...
            'ArtistCount',artists,'LeafCount',leaves);
        function execute()
            result = operation();
        end
    end
end

function result = roundtrip(file,payload)
mdemat.tools.io.writeRecord(file,struct('Payload',{payload}),overwrite=true);
value = mdemat.tools.io.readRecord(file);
description = dir(file);
result = struct('Value',{value.Payload},'SerializedBytes',description.bytes, ...
    'LeafCount',leafCount(payload));
end

function result = constructCollection(view)
result = struct();
for index = 1:8
    result.(sprintf('D%02d',index)) = struct('Type','numeric', ...
        'Data',struct('x',view,'v',view,'u2',view));
end
end

function result = updateCollection(original)
result = original;
result.D01.Data.u2.Values(2,1) = result.D01.Data.u2.Values(2,1)+1;
end

function count = leafCount(value)
if iscell(value)
    count = sum(cellfun(@leafCount,value),'all');
elseif isstruct(value)
    count = sum(cellfun(@leafCount,struct2cell(value)),'all');
else
    count = 1;
end
end

function result = solveOnce(A,b,recordHash)
solver = mdemat.common.linalg.LinearSolver(nnzTh=Inf,recordHash=recordHash);
value = solver.solve(A=A,b=b);
result = struct('Value',value,'HashRecorded',~isempty(solver.LhsHash));
end

function result = timerIntervals()
timer = mdemat.tools.analyze.Timer(verbose=false);
names = {'assemble','solve','assess'};
for index = 1:30
    timer.start(record=names{mod(index-1,3)+1});
    timer.stop();
end
result = timer.snapshot();
end

function valid = validTimer(result)
valid = ~result.IsActive && numel(result.Records)==3 && ...
    all([result.Records.Duration]>=0) && result.Total==sum([result.Records.Duration]);
end

function result = drawGraph(axis,vertices,edges,batch)
handles = gobjects(0,1);
try
    if batch
        handles = mdemat.data.mesh.internal.drawEdges(axis,vertices,edges,'k',1);
    else
        handles = gobjects(size(edges,1),1);
        for index = 1:size(edges,1)
            points = vertices(edges(index,:),:);
            handles(index) = line(axis,points(:,1),points(:,2),'Color','k','LineWidth',1);
        end
    end
    drawnow;
    if batch
        x = reshape(handles.XData,3,[]); y = reshape(handles.YData,3,[]);
        endpoints = [x(1,:).',y(1,:).',x(2,:).',y(2,:).'];
    else
        endpoints = zeros(size(edges,1),4);
        for index = 1:numel(handles)
            endpoints(index,:) = [handles(index).XData(1),handles(index).YData(1), ...
                handles(index).XData(2),handles(index).YData(2)];
        end
    end
    result = struct('Endpoints',endpoints,'ArtistCount',numel(handles));
catch exception
    delete(handles(isgraphics(handles)));
    rethrow(exception);
end
delete(handles(isgraphics(handles)));
end

function removeFigure(handle,previous)
if isgraphics(handle), delete(handle); end
if isempty(previous) || isgraphics(previous,'figure'), set(groot,'CurrentFigure',previous); end
end

function removeFolder(folder)
if ~isfolder(folder), return; end
[ok,message] = rmdir(folder,'s');
if ~ok, warning('mdemat:benchmark:toolsStudy:CleanupFailed','%s',message); end
end

function result = description(value)
stored = whos('value');
result = struct('Class',class(value),'Size',size(value),'Sparse',issparse(value), ...
    'Complex',~isreal(value),'Nonzeros',nnz(value),'PayloadBytes',stored.bytes);
end
