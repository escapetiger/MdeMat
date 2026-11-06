function config = studyConfig(config, options)
%STUDYCONFIG Explicit observers and a bounded version of a retained study.
%   Only research entrypoints call this helper. Numerical problem/solve APIs
%   continue to construct no observers. Smoke changes resolution and horizon,
%   while retaining the original physical callbacks and numerical method.
if ~options.Smoke, return; end
config.nx = repmat(4, 1, config.nDims);
config.xBasisOrder = 2;
config.verbose = 0;
% Steady studies have no time configuration. Keep their input schema exact.
if isfield(config, 'tFinal') || isfield(config, 'tOdeIntName')
    config.tFinal = 0.002;
    config.dt = 0.001;
    config.cfl = [];
    config.ckptTimeStamps = [];
end
if isfield(config, 'nh'), config.nh = min(config.nh, 4); end
if isfield(config, 'E')
    if isfield(config, 'nu'), config.nu = min(config.nu, 2); else, config.nu = 1; end
    if isfield(config, 'nv'), config.nv = min(config.nv, 4); else, config.nv = 4; end
end
end
