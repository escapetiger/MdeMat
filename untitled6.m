clear; clc; close all;

%% 1) Nodes on the unit circle
N = 20;
theta_nodes = linspace(0,2*pi,N+1)'; theta_nodes(end) = [];
x_nodes = [cos(theta_nodes), sin(theta_nodes)];

% Target function: smooth trend + high-frequency oscillation
f = @(theta) cos(theta) + 0.5*sin(10*theta);
f_nodes = f(theta_nodes);

%% 2) Fourier basis (smooth trend)
M = 3; % up to 3rd harmonic
Phi_nodes = zeros(N,2*M+1);
for k = 1:M
    Phi_nodes(:,k) = cos(k*theta_nodes);
    Phi_nodes(:,M+k) = sin(k*theta_nodes);
end
Phi_nodes(:,end) = ones(N,1); % constant term
M_total = size(Phi_nodes,2);

%% 3) Nodal basis (identity = Lagrange)
H = eye(N);

%% 4) Orthogonal system
A_ortho = [Phi_nodes H; Phi_nodes' zeros(M_total,M_total)];
rhs_ortho = [f_nodes; zeros(M_total,1)];
coef_ortho = A_ortho\rhs_ortho;
a_ortho = coef_ortho(1:M_total);
b_ortho = coef_ortho(M_total+1:end);

%% 5) Non-orthogonal system
A_no = [Phi_nodes H];
coef_no = A_no\f_nodes;
a_no = coef_no(1:M_total);
b_no = coef_no(M_total+1:end);

%% 6) Evaluation
theta_eval = linspace(0,2*pi,500).';
Phi_eval = zeros(length(theta_eval),M_total);
for k = 1:M
    Phi_eval(:,k) = cos(k*theta_eval);
    Phi_eval(:,M+k) = sin(k*theta_eval);
end
Phi_eval(:,end) = ones(length(theta_eval),1);

H_eval = zeros(length(theta_eval),N);
for j = 1:N
    ej = zeros(N,1); ej(j)=1;
    H_eval(:,j) = lagrange_basis_1D(theta_nodes, ej, theta_eval);
end

f_ortho = Phi_eval*a_ortho + H_eval*b_ortho;
g_ortho = Phi_eval*a_ortho; % smooth trend
h_ortho = H_eval*b_ortho;   % local residual

f_no = Phi_eval*a_no + H_eval*b_no;
g_no = Phi_eval*a_no;
h_no = H_eval*b_no;

f_true = f(theta_eval);

%% 7) Plot
figure('Name','Orthogonal vs Non-Orthogonal Hybrid');
plot(theta_eval,f_true,'k-','LineWidth',1.5); hold on;
plot(theta_eval,f_ortho,'r--','LineWidth',1.4);
plot(theta_eval,f_no,'b-.','LineWidth',1.4);
legend('f true','Orthogonal','Non-orthogonal'); grid on;
xlabel('\theta'); ylabel('f(\theta)'); title('Hybrid Approximation Comparison');

figure('Name','Contributions: Smooth vs Local (Orthogonal)');
plot(theta_eval,g_ortho,'r--','LineWidth',1.4); hold on;
plot(theta_eval,h_ortho,'r:','LineWidth',1.4);
legend('Smooth (g)','Local (h)'); xlabel('\theta'); ylabel('Contribution'); grid on;

figure('Name','Contributions: Smooth vs Local (Non-Orthogonal)');
plot(theta_eval,g_no,'b--','LineWidth',1.4); hold on;
plot(theta_eval,h_no,'b:','LineWidth',1.4);
legend('Smooth (g)','Local (h)'); xlabel('\theta'); ylabel('Contribution'); grid on;

%% -----------------------------
% Helper: barycentric Lagrange 1D
function y = lagrange_basis_1D(x_nodes, values, x_eval)
    N = length(x_nodes);
    w = ones(N,1);
    for j=1:N
        for k=[1:j-1,j+1:N]
            w(j) = w(j)/(sin((x_nodes(j)-x_nodes(k))/2));
        end
    end
    num = zeros(size(x_eval));
    den = zeros(size(x_eval));
    for j=1:N
        diff = sin((x_eval - x_nodes(j))/2);
        term = w(j)./diff;
        num = num + values(j).*term;
        den = den + term;
    end
    y = num./den;
    % handle exact nodes
    for j=1:N
        mask = abs(mod(x_eval - x_nodes(j)+pi,2*pi)-pi)<1e-14;
        y(mask) = values(j);
    end
end
