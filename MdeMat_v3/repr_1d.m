%% Hybrid approximation comparison (orthogonalized vs not)
clear; clc;

% -------------------------------
% 1) Target function and nodes
% -------------------------------
f = @(x) 1./(1+25*x.^2) + 0.02*cos(18*x);  % test function
N = 10;                                   % number of nodes for nodal part
x_nodes = lgwt(N,-1,1);                   % Gauss–Legendre nodes (column vector)
f_nodes = f(x_nodes(:));                  % column vector (N x 1)

% Legendre polynomial degree
M = 1;  % number of Legendre terms

% -------------------------------
% 2) Build Legendre matrix at nodes
% -------------------------------
P = zeros(N,M);
for k = 0:M-1
    P(:,k+1) = legendreP(k,x_nodes);
end

% Evaluation grid
x_eval = linspace(-1,1,800).';
P_eval = zeros(length(x_eval),M);
for k = 0:M-1
    P_eval(:,k+1) = legendreP(k,x_eval);
end
f_true = f(x_eval).';

%% -----------------------------
% 3) Orthogonalized Hybrids
%% -----------------------------
% ---- Legendre + Lagrange (orthogonalized) ----
H_L = eye(N);
A_L_ortho = [P H_L; P' zeros(M,M)];
rhs_L_ortho = [f_nodes; zeros(M,1)];
sol_L_ortho = A_L_ortho\rhs_L_ortho;
a_L_ortho = sol_L_ortho(1:M);
b_L_ortho = sol_L_ortho(M+1:end);

H_eval_L = zeros(length(x_eval),N);
for j = 1:N
    ej = zeros(N,1); ej(j)=1;
    H_eval_L(:,j) = lagrange_basis(x_nodes, ej, x_eval);
end
f_L_ortho = P_eval*a_L_ortho + H_eval_L*b_L_ortho;
g_L_ortho = P_eval*a_L_ortho;
h_L_ortho = H_eval_L*b_L_ortho;

% ---- Legendre + RBF (orthogonalized) ----
epsilon = 2; 
rbf = @(r) exp(-(epsilon*r).^2);
D = abs(x_nodes - x_nodes.'); % N x N
H_R = rbf(D);

A_R_ortho = [P H_R; P' zeros(M,M)];
rhs_R_ortho = [f_nodes; zeros(M,1)];
sol_R_ortho = A_R_ortho\rhs_R_ortho;
a_R_ortho = sol_R_ortho(1:M);
b_R_ortho = sol_R_ortho(M+1:end);

D_eval = abs(x_eval(:) - x_nodes(:).'); % nEval x N
H_eval_R = rbf(D_eval);

f_R_ortho = P_eval*a_R_ortho + H_eval_R*b_R_ortho;
g_R_ortho = P_eval*a_R_ortho;
h_R_ortho = H_eval_R*b_R_ortho;

%% -----------------------------
% 4) Non-Orthogonal Hybrids
%% -----------------------------
% ---- Legendre + Lagrange (no orthogonality) ----
A_L_no = [P H_L];
sol_L_no = A_L_no\f_nodes;
a_L_no = sol_L_no(1:M);
b_L_no = sol_L_no(M+1:end);
f_L_no = P_eval*a_L_no + H_eval_L*b_L_no;

% ---- Legendre + RBF (no orthogonality) ----
A_R_no = [P H_R];
sol_R_no = A_R_no\f_nodes;
a_R_no = sol_R_no(1:M);
b_R_no = sol_R_no(M+1:end);
f_R_no = P_eval*a_R_no + H_eval_R*b_R_no;

%% -----------------------------
% 5) Error analysis
%% -----------------------------
err_L_ortho = abs(f_true - f_L_ortho);
err_R_ortho = abs(f_true - f_R_ortho);
err_L_no = abs(f_true - f_L_no);
err_R_no = abs(f_true - f_R_no);

fprintf('Max errors:\n');
fprintf('Legendre+Lagrange: ortho %.3e, no ortho %.3e\n', max(err_L_ortho), max(err_L_no));
fprintf('Legendre+RBF     : ortho %.3e, no ortho %.3e\n', max(err_R_ortho), max(err_R_no));

%% -----------------------------
% 6) Plots
%% -----------------------------
% Main comparison
figure('Name','Hybrid Approximations');
plot(x_eval,f_true,'k-','LineWidth',1.6); hold on;
plot(x_eval,f_L_ortho,'r--','LineWidth',1.4);
plot(x_eval,f_R_ortho,'b-.','LineWidth',1.4);
plot(x_eval,f_L_no,'r:','LineWidth',1.2);
plot(x_eval,f_R_no,'b:','LineWidth',1.2);
plot(x_nodes,f_nodes,'ko','MarkerFaceColor','k');
legend('f true','Legendre+Lagrange ortho','Legendre+RBF ortho',...
    'Legendre+Lagrange no ortho','Legendre+RBF no ortho','nodes');
xlabel('x'); ylabel('f(x)'); grid on; title('Hybrid Approximation Comparison');

% Absolute errors
figure('Name','Absolute Errors (log scale)');
semilogy(x_eval,err_L_ortho,'r--','LineWidth',1.4); hold on;
semilogy(x_eval,err_R_ortho,'b-.','LineWidth',1.4);
semilogy(x_eval,err_L_no,'r:','LineWidth',1.2);
semilogy(x_eval,err_R_no,'b:','LineWidth',1.2);
legend('Legendre+Lagrange ortho','Legendre+RBF ortho','Legendre+Lagrange no ortho','Legendre+RBF no ortho');
xlabel('x'); ylabel('|f - f_{approx}|'); grid on; title('Absolute Errors (log scale)');

%% -----------------------------
% 7) Helper: barycentric Lagrange
%% -----------------------------
function y = lagrange_basis(x_nodes, values, x_eval)
    N = length(x_nodes);
    w = ones(N,1);
    for j=1:N
        for k=[1:j-1,j+1:N]
            w(j) = w(j)/(x_nodes(j)-x_nodes(k));
        end
    end
    num = zeros(size(x_eval));
    den = zeros(size(x_eval));
    for j=1:N
        diff = x_eval - x_nodes(j);
        term = w(j)./diff;
        num = num + values(j).*term;
        den = den + term;
    end
    y = num./den;
    % handle exact nodes
    for j=1:N
        mask = abs(x_eval - x_nodes(j))<1e-14;
        y(mask) = values(j);
    end
end
