%==========================================================================
% FileName: ex01_interpolation_basis_compilation
% Author: Yi Cai
% Contact: yicaim@stu.xmu.edu.cn
% Date: 2025-05-27
% Description:
%   Compile interpolation basis functions.
%==========================================================================

clear; clc;

%% 1. Lagrange Basis with Gauss-Legendre nodes on [-1/2, 1/2]

n = 5;

fprintf('Compiling Lagrange basis on Gauss-Legendre nodes...\n');
GL = approx.integrate.GaussLegendreRule();
for i = 1:n
    compiler = core.symbolic.InterpolationPolynomialBasisCompiler();
    [nodes, ~] = GL.generate(i, false, -0.5, 0.5);
    compiler.setFunctions(nodes, 'lagrange');
    compiler.compile(sprintf('lagrange_unit_gauss_legendre_%d', i), 'lagrange', max(1, i-1));
end

%% 2. Lagrange Basis with Gauss-Lobatto nodes on [-1/2, 1/2]

n = 5;

fprintf('Compiling Lagrange basis on Gauss-Lobatto nodes...\n');
GL = approx.integrate.GaussLobattoRule();
for i = 1:n
    compiler = core.symbolic.InterpolationPolynomialBasisCompiler();
    [nodes, ~] = GL.generate(i, false, -0.5, 0.5);
    compiler.setFunctions(nodes, 'lagrange');
    compiler.compile(sprintf('lagrange_unit_gauss_lobatto_%d', i), 'lagrange', max(1, i-1));
end
