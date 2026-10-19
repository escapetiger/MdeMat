clear; clc; close all;

%% 1) Nodes
N = 16;                     % number of nodes
x_nodes = linspace(-1,1,N)';  % equidistant nodes

%% 2) Test function
f = @(x) 1./(1 + 25*x.^2);    % Runge function

f_nodes = f(x_nodes);

%% 3) Lagrange basis interpolation
x_eval = linspace(-1,1,500)';
L_poly = zeros(N,length(x_eval));
for j = 1:N
    % Lagrange polynomial L_j(x)
    Lj = ones(size(x_eval));
    for k = 1:N
        if k ~= j
            Lj = Lj .* (x_eval - x_nodes(k)) / (x_nodes(j) - x_nodes(k));
        end
    end
    L_poly(j,:) = Lj;
end
f_lagrange = f_nodes' * L_poly;

%% 4) Gaussian RBF
epsilon = 5; % shape parameter
H = exp(-(pdist2(x_nodes,x_nodes)*epsilon).^2);

% Construct cardinal RBFs
C = H \ eye(N); % weights for cardinal RBFs
rbf_eval = exp(-(pdist2(x_eval,x_nodes)*epsilon).^2);
f_rbf = rbf_eval * C * f_nodes;

%% 5) Plot comparison
figure('Name','Interpolation Comparison');
plot(x_eval,f(x_eval),'k-','LineWidth',1.5); hold on;
plot(x_eval,f_lagrange,'r--','LineWidth',1.4);
plot(x_eval,f_rbf,'b-.','LineWidth',1.4);
plot(x_nodes,f_nodes,'ko','MarkerFaceColor','k');
legend('True function','Lagrange','Cardinal RBF','Nodes');
xlabel('x'); ylabel('f(x)'); grid on;
title('Lagrange vs Cardinal Gaussian RBF Interpolation');

%% 6) Visualize basis functions
figure('Name','Basis Functions Comparison');
for j = 1:N
    plot(x_eval,L_poly(j,:),'r--','LineWidth',1); hold on;
    plot(x_eval,rbf_eval*C(:,j),'b-.','LineWidth',1);
end
legend('Lagrange','Cardinal RBF');
xlabel('x'); ylabel('Basis function value'); grid on;
title('Lagrange vs Cardinal Gaussian RBF Basis Functions');
