%% Configuration and foundational tools
% Config is a value map. Every declared key owns a Value/Validator pair;
% updates must be assigned back. The predicate is read-only and returns one
% logical scalar. No validator is represented by [].
%
%   cfg = mdemat.tools.environment.Config();
%   validStep = @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0);
%   cfg = cfg.define("dt", [], validStep);
%   cfg = cfg.set("dt", 0.01);
%   values = cfg.toStruct();
%
% Definitions can contain deferred defaults. validate(keys) checks the
% requirements of a particular operation; it does not interpret [] as a
% universal missing value. update validates an entire scalar struct patch
% before returning a changed Config. Copying Config isolates entry changes;
% embedded handles retain normal MATLAB shallow sharing.
%
%   [patch, info] = mdemat.tools.environment.readConfig(configFile);
%   cfg = cfg.update(patch); % every patch key must already be declared
%
% readConfig does not execute MATLAB, register predicates or change path/cwd.
% It accepts OPTION=name=literal, quoted text, scalar arithmetic from a small
% mathematical allowlist, and rectangular arrays/cells. Use parentheses for
% arithmetic array entries containing spaces. Duplicate OPTION keys fail.
% MATLABPATH/SUBDIRS are source-relative descriptions in info.Paths, separate
% from values. info.File and info.Lines provide source diagnostics.
%
%% Scheme ownership
% Schemes compose Config. Known inputs use the model owner's predicate even
% when a caller supplies a different predicate in an explicit Config. A
% scalar struct supplies values for registered keys; an explicit Config may
% declare additional extension keys. Derived/resource names remain reserved.
%
%   cfg = mdemat.tools.environment.Config();
%   cfg = cfg.define("tFinal", 1);
%   cfg = cfg.define("experimentId", "trialA");
%   scheme = mdemat.internal.execution.OdeScheme(config=cfg);
%
% Construction freezes numerical bindings alongside the selected method.
% Later numerical changes require a new Scheme. Partial instances may keep
% deferred defaults for operations such as explicit restart; initialize must
% validate its own required keys. Config is a read-only plain value snapshot.
% Validators, model callbacks and observer handles are not persistence data.
%
% setObserver/getObserver control timer and visualizer resources separately.
% setObservation("plot", patch) and setObservation("analysis", patch) shallow
% merge explicit observation settings without changing numerical bindings.
% PointsPerDimension is a per-axis count; scalar zero disables sampling.
%
%% Scoped wall-clock timing
% Timer is the only timing handle. Records/find/snapshot return ordinary
% values: edits and later intervals do not change an earlier snapshot.
%
%   timer = mdemat.tools.analyze.Timer(verbose=false);
%   [~, interval] = timer.start(record="solve");
%   guard = onCleanup(@() timer.cancel(interval));
%   % Perform the operation here.
%   timer.stop(interval=interval);
%   clear guard
%   timing = timer.snapshot();
%
% Duration includes stopped intervals only. cancel discards its matching
% active interval; stale cleanup cannot stop newer work. Named durations
% accumulate, intervals are mutually exclusive, and report prints wall time
% in seconds. Reset clears durations and cancels active work with a warning.
%
%% Numerical identities and environment
% tools.base.hash uses the class, shape, sparse/complex flags and little-endian
% numeric bytes. Sparse entries include uint64 row/column coordinates in find
% order; real bytes precede imaginary bytes. This is exact representation
% identity, including signed zeros and NaN payloads, not numerical closeness.
% hex is lowercase char text; uint8/double are rows of digest byte values.
% Hashing requires a JVM; LinearSolver leaves diagnostic hashing disabled by
% default (recordHash=true opts in). Equality controls factorization reuse.
%
% toolboxRoot locates installed resources independently of cwd. projectRoot
% without arguments accepts only this function's source checkout layout;
% an explicit directory searches its ancestors. versionInfo validates trusted
% installed metadata. describe excludes personal/license information and only
% lists installed products on explicit request; inventory is not dependency
% or license-availability evidence.
