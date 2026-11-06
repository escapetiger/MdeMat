%% EX06_SPHERICAL_HARMONIC_BASIS Generate spherical harmonic orthogonal basis functions.
%
%   This script generates spherical harmonic basis functions up to a
%   specified maximum degree. Configure parameters below and run the
%   script.

%% Configuration Parameters
maxDegree = 5; % Maximum degree l
variables = [sym('mu', 'real'), sym('phi', 'real')]; % mu=cos(theta), azimuth phi
isNormalized = true; % Orthonormal basis flag
symbolicDirectory = fileparts(mfilename('fullpath'));
fileName = fullfile(symbolicDirectory, 'metadata', 'normal_spherical_harmonic');
maxDerivOrder = 0; % Maximum derivative order for handles

%% Generate Spherical Harmonic Basis Functions
fprintf('=== Spherical Harmonic Basis Generation ===\n');
fprintf('Generating spherical harmonic basis up to degree %d...\n', ...
    maxDegree);

basis = mdemat.common.symbolic.SphericalHarmonicBasis(variables, maxDegree, ...
    isNormalized=isNormalized);
basis.compile(maxDerivOrder);

%% Display Results
fprintf('Generated %d basis functions:\n', basis.NCodims);

for i = 1:basis.NCodims
    fprintf('Y_%d(mu,phi) = %s\n', i-1, char(basis.Expressions(i)));
end
fprintf('\n');

fprintf('Basis Properties:\n');
fprintf('  Max Degree: %d\n', basis.MaxDegree);
fprintf('  IsNormalized: %s\n', string(basis.IsNormalized));
fprintf('  Input Dims: %d\n', basis.NDims);
fprintf('  Output Dims: %d\n', basis.NCodims);
fprintf('\n');

%% Save Results (if specified)
if ~isempty(fileName)
    basis.save(fileName);
    fprintf('Basis saved to: %s\n', [fileName, basis.MetaExt]);
end
