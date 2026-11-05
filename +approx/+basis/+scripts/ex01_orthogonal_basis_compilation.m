%==========================================================================
% FileName: ex01_orthogonal_basis_compilation
% Author: Yi Cai
% Contact: yicaim@stu.xmu.edu.cn
% Date: 2025-05-27
% Description:
%   Compile orthogonal basis functions.
%==========================================================================

clear; clc;

%% 1. Legendre Basis on [-1/2, 1/2]

n = 4;
for i = 1:n
    compiler = core.symbolic.OrthogonalPolynomialBasisCompiler();
    fprintf('Compiling monic Legendre basis on [-1/2, 1/2]...\n');
    compiler.setFunctions(i-1, [-0.5, 0.5], 'monic_legendre');
    compiler.compile(sprintf('monic_legendre_unit_%d', i), 'monic_legendre', max(1, i-1));
end

%% 2. Legendre Basis on [-1, 1]

n = 16;

for i = 1:n
    compiler = core.symbolic.OrthogonalPolynomialBasisCompiler();
    fprintf('Compiling Legendre basis on [-1, 1]...\n');
    compiler.setFunctions(i-1, [-1, 1], 'legendre');
    compiler.compile(sprintf('legendre_canonical_%d', i), 'legendre', max(1, i-1));
end


%% 3. Fourier basis on [0, 2*pi]

n = 16;

for i = 1:n
    compiler = core.symbolic.FourierBasisCompiler();
    fprintf('Compiling real-valued Fourier basis on [0, 2*pi]...\n');
    compiler.setFunctions(i, [0, 2*pi], 'fourier');
    compiler.compile(sprintf('fourier_canonical_%d', i), 'fourier', 0);
end


%% 3. Spherical harmonic basis on [-1, 1]x[0, 2*pi]

n = 16;

for i = 1:n
    compiler = core.symbolic.SphericalHarmonicBasisCompiler({'mu', 'theta'});
    fprintf('Compiling real-valued Spherical harmonic basis on [-1, 1]x[0, 2*pi]...\n');
    compiler.setFunctions(i, 'spherical_harmonic');
    compiler.compile(sprintf('spherical_harmonic_%d', i), 'spherical_harmonic', 0);
end