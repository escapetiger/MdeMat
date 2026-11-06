function StateFiles
%STATEFILES Versioned state files and completed-step restart.
%   State.save(FILE) stores numeric DOFs, coefficients and observation history
%   with explicit space descriptions. It refuses overwrite and missing parent
%   folders; opt in using overwrite=true and createDirectory=true. Both MAT and
%   JSON use mdemat.tools.io's fixed data schema. Old opaque MATLAB class saves
%   are refused rather than automatically deserialized.
%
%   state.save(file, createDirectory=true)
%   restored = mdemat.common.fem.SpatialState.load(file)
%   State.load reads State snapshots only. To extract a State for observation
%   from an execution checkpoint, select both operations explicitly:
%   checkpoint = mdemat.internal.execution.readCheckpoint(file);
%   restored = mdemat.common.fem.SpatialState.fromRecord(checkpoint.state);
%
%   Every built-in state has a matching class-specific load/fromRecord method.
%   Basis order, normalization, complex data, physical quadrature weights,
%   periodic mesh maps, velocity sum/scaling and macro-micro decomposition are
%   preserved. Solver/mesh/SL/visualization caches are rebuilt. Compiled basis
%   definitions are checked against trusted installed resources; a changed
%   definition raises ResourceChanged rather than reinterpreting old DOFs.
%
%   Coefficient callbacks require explicit omission at save and binding at load:
%   state.save(file, omitCallbacks=true)
%   restored = mdemat.common.fem.SpatialState.load(file, callbacks=struct('f', @sin))
%   Supply exactly the callback names named in the saved record. Callback
%   meaning and external data are the caller's responsibility.
%
%   scheme.saveCheckpoint(file, state) is called at an accepted step after
%   postStep/advance. It also stores BDF solution/mass/step histories and every
%   split factor's workspace. Model callbacks remain caller inputs:
%   state = scheme.loadCkptState(file, explicitStatePrototype)
%   state = scheme.run(state, resume=true)
%   The explicit configured scheme/model must match; loading validates new
%   objects before committing its ODE workspace. resume=true rebuilds equation
%   operators while retaining the accepted state/history. The automatic writer
%   uses ckptOverwrite=false unless explicitly configured otherwise.
%   Pure observation options do not change model identity. Saved v1 model
%   descriptions are projected before comparison; files are never rewritten.
%
%   See also mdemat.tools.io.writeRecord, mdemat.tools.io.readRecord,
%   mdemat.common.fem.describeSpace, mdemat.common.fem.restoreSpace,
%   mdemat.common.ode.snapshotIntegrator, mdemat.common.ode.restoreIntegrator.
end
