clear; clc; close all;

%% 1) Nodes on the sphere
% N = 50;
% theta = acos(2*rand(N,1)-1);    % polar angle
% phi   = 2*pi*rand(N,1);         % azimuth
% x_nodes = [sin(theta).*cos(phi), sin(theta).*sin(phi), cos(theta)]; % N x 3

N = 64;
i = (0:N-1)';
phi = 2*pi*i/((1+sqrt(5))/2);      % golden angle
theta = acos(1 - 2*(i+0.5)/N);     % latitude
x_nodes = [sin(theta).*cos(phi), sin(theta).*sin(phi), cos(theta)];

% n1 = 8;
% t1 = linspace(0, 1, n1).';
% n2 = 8;
% t2 = linspace(0, 1-1/n2, n2).';
% theta = acos(2*t1-1);    % polar angle
% phi   = 2*pi*t2;         % azimuth
% [theta, phi] = ndgrid(theta, phi);
% theta = theta(:);
% phi = phi(:);
% x_nodes = [sin(theta).*cos(phi), sin(theta).*sin(phi), cos(theta)]; % N x 3
% [x_nodes, I, ~] = unique(x_nodes, 'rows', 'stable');
% theta = theta(I);
% phi = phi(I);
% N = size(x_nodes, 1);

% Target function
f = @(x) x(:,1).^2 - x(:,2).^2 + 0.1*x(:,3);
f_nodes = f(x_nodes);

%% 2) Spherical harmonics (real) up to degree L
L = 3;
M = (L+1)^2;
Y_nodes = zeros(N,M);
k = 1;
for l = 0:L
    for m = -l:l
        Y_nodes(:,k) = real_sph_harm(l,m,theta,phi);
        k = k+1;
    end
end

%% 3) RBF part
epsilon = 3;
% rbf = @(r) exp(-(epsilon*r).^2);
% D = pdist2(x_nodes, x_nodes);
% H_R0 = rbf(D);

rbf = approx.basis.RbfKernelBasis(x_nodes.', 'gaussian', epsilon);
H_R = rbf.evaluate(x_nodes.').';

%% 4) Lagrange basis on sphere (identity for nodes)
H_Lagr = eye(N);

%% 5) Augmented system: SPH + RBF
A_RBF = [Y_nodes H_R; Y_nodes' zeros(M,M)];
rhs_RBF = [f_nodes; zeros(M,1)];
coef_RBF = A_RBF\rhs_RBF;
a_RBF = coef_RBF(1:M);
b_RBF = coef_RBF(M+1:end);

%% 6) Augmented system: SPH + Lagrange
A_Lagr = [Y_nodes H_Lagr; Y_nodes' zeros(M,M)];
rhs_Lagr = [f_nodes; zeros(M,1)];
coef_Lagr = A_Lagr\rhs_Lagr;
a_Lagr = coef_Lagr(1:M);
b_Lagr = coef_Lagr(M+1:end);

%% 7) Evaluation grid
[theta_grid, phi_grid] = meshgrid(linspace(0,pi,50), linspace(0,2*pi,100));
x_eval = [sin(theta_grid(:)).*cos(phi_grid(:)), ...
          sin(theta_grid(:)).*sin(phi_grid(:)), ...
          cos(theta_grid(:))];

% Evaluate spherical harmonics
Y_eval = zeros(size(x_eval,1),M);
k = 1;
for l = 0:L
    for m = -l:l
        Y_eval(:,k) = real_sph_harm(l,m,acos(x_eval(:,3)),atan2(x_eval(:,2),x_eval(:,1)));
        k = k+1;
    end
end

% RBF evaluation
% rbf = @(r) exp(-(epsilon*r).^2);
% D_eval = pdist2(x_eval, x_nodes);
% H_eval_RBF = rbf(D_eval);

H_eval_RBF = rbf.evaluate(x_eval.').';

% Lagrange evaluation (barycentric)
H_eval_Lagr = zeros(size(x_eval,1),N);
for j = 1:N
    ej = zeros(N,1); ej(j)=1;
    H_eval_Lagr(:,j) = lagrange_basis_sphere(x_nodes, ej, x_eval);
end

%% 8) Approximation
f_RBF = Y_eval*a_RBF + H_eval_RBF*b_RBF;
f_Lagr = Y_eval*a_Lagr + H_eval_Lagr*b_Lagr;
f_true = f(x_eval);

%% 9) Max errors
err_RBF = abs(f_true - f_RBF);
err_Lagr = abs(f_true - f_Lagr);

fprintf('Max error SPH + RBF: %.3e\n', max(err_RBF));
fprintf('Max error SPH + Lagr: %.3e\n', max(err_Lagr));

%% 10) Plots
f_RBF_grid = reshape(f_RBF,size(theta_grid));
f_Lagr_grid = reshape(f_Lagr,size(theta_grid));
f_true_grid = reshape(f_true,size(theta_grid));

figure;
subplot(1,3,1);
surf(sin(theta_grid).*cos(phi_grid), sin(theta_grid).*sin(phi_grid), cos(theta_grid), f_true_grid);
title('True function'); shading interp; axis equal; colorbar;

subplot(1,3,2);
surf(sin(theta_grid).*cos(phi_grid), sin(theta_grid).*sin(phi_grid), cos(theta_grid), f_RBF_grid);
title('SPH + RBF'); shading interp; axis equal; colorbar;

subplot(1,3,3);
surf(sin(theta_grid).*cos(phi_grid), sin(theta_grid).*sin(phi_grid), cos(theta_grid), f_Lagr_grid);
title('SPH + Lagrange'); shading interp; axis equal; colorbar;

%% -----------------------------
% Helper: real spherical harmonics
function Y = real_sph_harm(l,m,theta,phi)
    P = legendre(l, cos(theta)); 
    if m > 0
        K = sqrt(2) * P(m+1,:)' .* cos(m*phi);
    elseif m < 0
        K = sqrt(2) * P(abs(m)+1,:)' .* sin(abs(m)*phi);
    else
        K = P(1,:)';
    end
    Y = K;
end

%% -----------------------------
% Helper: barycentric Lagrange on sphere
function y = lagrange_basis_sphere(nodes, values, x_eval)
    % nodes: N x 3, x_eval: nEval x 3
    N = size(nodes,1);
    w = ones(N,1);
    for j=1:N
        for k=[1:j-1,j+1:N]
            w(j) = w(j)/(norm(nodes(j,:)-nodes(k,:)));
        end
    end
    num = zeros(size(x_eval,1),1);
    den = zeros(size(x_eval,1),1);
    for j=1:N
        diff = sqrt(sum((x_eval - nodes(j,:)).^2,2));
        term = w(j)./diff;
        num = num + values(j)*term;
        den = den + term;
    end
    y = num./den;
    % handle exact node matches
    for j=1:N
        mask = sqrt(sum((x_eval - nodes(j,:)).^2,2))<1e-14;
        y(mask) = values(j);
    end
end
