function pairs = gaussianSweep(name)
%GAUSSIANSWEEP Select a parameter sweep for the Gaussian radiation example.
%   PAIRS = gaussianSweep() returns the original varying-nv sweep.
%   PAIRS = gaussianSweep("varyNv") holds nu=6 and varies nv.
%   PAIRS = gaussianSweep("varyNu") holds nv=4 and varies nu.
%
%   Each row contains [nu, nv]. The result is a 4-by-2 double array.
%   The two presets preserve the original example and its desktop variant.
%
%   See also ex03_gaussian.

arguments
    name (1, 1) string = "varyNv"
end

switch name
    case "varyNv"
        pairs = [6, 0; 6, 4; 6, 8; 6, 16];
    case "varyNu"
        pairs = [1, 4; 4, 4; 8, 4; 16, 4];
    otherwise
        error('physics:radiation:InvalidGaussianSweep', ...
            'Unknown Gaussian sweep "%s". Use "varyNv" or "varyNu".', name);
end
end
