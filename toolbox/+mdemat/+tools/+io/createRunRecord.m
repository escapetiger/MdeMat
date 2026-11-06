function record = createRunRecord(solution, solveReport, options)
%CREATERUNRECORD Keep solution, execution, accuracy and provenance distinct.
%   R = createRunRecord(SOLUTION, SOLVEREPORT) accepts plain value structs.
%   An omitted validation report is explicitly UNVERIFIED; successful execution
%   does not imply accuracy. A supplied validation report declares Status,
%   Criterion, Measurements, Reference, Bounds and Reason. PASS/FAIL require
%   nonempty criterion/reference/bounds/measurements supplied by the assessor.
%
%   inputs, timings, utc (yyyy-MM-ddTHH:mm:ssZ), seed and codeVersion are
%   explicit caller inputs. Empty provenance fields mean not supplied. This
%   function does not inspect Git, read clocks/RNG or run/assess a solver.
%   environment defaults to the native MATLAB version/architecture inventory;
%   installed products remain opt-in and never imply dependencies/licenses.
arguments
    solution (1, 1) struct
    solveReport (1, 1) struct
    options.validation (1, 1) struct = struct()
    options.inputs (1, 1) struct = struct()
    options.timings (1, 1) struct = struct()
    options.environment (1, 1) struct = struct()
    options.utc {mustBeTextScalar} = ''
    options.seed = []
    options.codeVersion {mustBeTextScalar} = ''
    options.limits (1,1) struct = struct()
end
identity = mdemat.tools.environment.versionInfo();
if identity.RunRecordSchemaVersion ~= 1 || identity.NumericRecordSchemaVersion ~= 1
    error('mdemat:io:UnsupportedSchema','This factory implements run/numeric record schema v1.');
end
validation = options.validation;
if isempty(fieldnames(validation))
    validation = struct('Status', 'UNVERIFIED', 'Criterion', '', ...
        'Measurements', struct(), 'Reference', '', 'Bounds', struct(), ...
        'Reason', 'No independent accuracy assessment was supplied.');
end
required = {'Status', 'Criterion', 'Measurements', 'Reference', 'Bounds', 'Reason'};
if ~isequal(sort(fieldnames(validation)), sort(required(:)))
    error('mdemat:io:InvalidValidationReport', 'Use the explicit validation report schema.');
end
for field = ["Status", "Criterion", "Reference", "Reason"]
    value = validation.(field);
    if ~((ischar(value) && (isrow(value) || isempty(value))) || ...
            (isstring(value) && isscalar(value) && ~ismissing(value)))
        error('mdemat:io:InvalidValidationReport', 'Validation text must be a scalar string/character row.');
    end
    validation.(field) = char(value);
end
if ~ismember(validation.Status, {'PASS', 'FAIL', 'UNVERIFIED'}) || ...
        ~isstruct(validation.Measurements) || ...
        ~isscalar(validation.Measurements) || ~isstruct(validation.Bounds) || ~isscalar(validation.Bounds)
    error('mdemat:io:InvalidValidationReport', 'Use the explicit validation report schema.');
end
if ~strcmp(validation.Status, 'UNVERIFIED') && (isempty(validation.Criterion) || ...
        isempty(validation.Reference) || isempty(fieldnames(validation.Measurements)) || ...
        isempty(fieldnames(validation.Bounds)))
    error('mdemat:io:InvalidValidationReport', ...
        'PASS/FAIL requires criterion, reference, measurements and bounds.');
end
utc = char(options.utc);
if ~isempty(utc)
    if isempty(regexp(utc, '^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$', 'once'))
        error('mdemat:io:InvalidProvenance', 'UTC must use yyyy-MM-ddTHH:mm:ssZ.');
    end
    try
        datetime(utc, InputFormat="yyyy-MM-dd'T'HH:mm:ss'Z'", TimeZone='UTC');
    catch exception
        failure = MException('mdemat:io:InvalidProvenance', 'Invalid UTC date/time.');
        throwAsCaller(addCause(failure, exception));
    end
end
seed = options.seed;
if ~isempty(seed) && (~isnumeric(seed) || ~isscalar(seed) || ~isreal(seed) || ...
        ~isfinite(seed) || seed < 0 || seed > 2^32-1 || seed ~= fix(seed))
    error('mdemat:io:InvalidProvenance', 'Seed must be an explicit integer in [0,2^32-1] or empty.');
end
environment = options.environment;
if isempty(fieldnames(environment)), environment = mdemat.tools.environment.describe(); end
record = struct('Format', 'mdemat.runRecord', 'SchemaVersion', identity.RunRecordSchemaVersion, ...
    'ToolboxVersion', identity.ToolboxVersion, ...
    'NumericRecordSchemaVersion', identity.NumericRecordSchemaVersion, ...
    'Inputs', options.inputs, 'Solution', solution, 'SolveReport', solveReport, ...
    'ValidationReport', validation, 'Timings', options.timings, 'Environment', environment, ...
    'Provenance', struct('UTC', utc, 'Seed', seed, 'CodeVersion', char(options.codeVersion)));
% Reuse the actual persistence grammar instead of inventing a second value validator.
packRecord(record,recordLimits(options.limits));
end
