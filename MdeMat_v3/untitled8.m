%% Full 1D x-v demo (corrected macro reconstruction)
% - Transport: partial_t f + v * partial_x f = (1/eps) L(f)
% - Collision: BGK-like linear L(f) = -nu * (f - P f)
% - Micro: cardinal Gaussian RBFs at velocity nodes (h_j = f_micro(v_j))
% - Projection P computed by quadrature on a separate fine grid

clear; close all; clc;

%% -------------------------
% Parameters
%% -------------------------
x_min = 0; x_max = 1;    Nx = 100;           % spatial grid
Nv    = 21;               % number of velocity nodes (RBF centers)
v_max = 6; v_min = -6;    % velocity domain for nodes
Nq    = 801;             % quadrature points in velocity
Tfinal = 0.5; dt = 5e-4; % time stepping
eps = 1e-2;              
nu  = 1/eps;             

% RBF parameters
epsilon_rbf = 0.8;       

% Macro basis (moments) phi_k(v) = [1, v, v^2]
K = 3;

%% -------------------------
% Grids
%% -------------------------
x = linspace(x_min, x_max, Nx)'; dx = x(2)-x(1);
v_nodes = linspace(v_min, v_max, Nv)';    % RBF centers

% quadrature grid for velocity integrals (finer, independent of RBF nodes)
vq = linspace(v_min*1.5, v_max*1.5, Nq)'; 
wq = ones(Nq,1); wq(2:end-1)=2; wq = wq*(vq(end)-vq(1))/(2*(Nq-1)); 

%% -------------------------
% Build cardinal RBFs
%% -------------------------
D_nn = abs(v_nodes - v_nodes.');          
Phi_nn = exp(-(epsilon_rbf * D_nn).^2);   

C = Phi_nn \ eye(Nv);                     

D_qn = abs(vq - v_nodes.');
Phi_qn = exp(-(epsilon_rbf * D_qn).^2);   
Psi_tilde_q = Phi_qn * C;                 

Id_check = (exp(-(epsilon_rbf * D_nn).^2) * C);
fprintf('||H*C - I||_inf = %.3e (should be ~0)\n', norm(Id_check - eye(Nv), inf));

%% -------------------------
% Macro basis and projection weights
%% -------------------------
Phi_macro_q = [ones(Nq,1), vq, vq.^2];    

M = (Phi_macro_q' .* wq') * Phi_macro_q;  
M_inv = inv(M);

W = (Psi_tilde_q' .* wq') * Phi_macro_q;  

Phi_macro_nodes = [ones(Nv,1), v_nodes, v_nodes.^2];   
P_node = Phi_macro_nodes * (M_inv * (W'));  

Lmat = -nu * (eye(Nv) - P_node);

%% -------------------------
% Initial condition
%% -------------------------
h0 = zeros(Nv,Nx);
for ix = 1:Nx
    h0(:,ix) = 0.2 * sin(2 * pi * v_nodes / v_max) .* exp(-((x(ix)-0.5)/0.15).^2);
end
h = h0;    % micro coefficients

%% -------------------------
% Pre-factor collision solver
%% -------------------------
Acol = eye(Nv) - (dt/eps) * Lmat;   
[Llu, Ulu, Plu] = lu(Acol);         

%% -------------------------
% Fine grid for plotting micro
%% -------------------------
vq_plot = linspace(v_min*1.2, v_max*1.2, 201)'; 
Phi_plot = exp(-(epsilon_rbf * abs(vq_plot - v_nodes')).^2) * C; 

%% -------------------------
% Time stepping
%% -------------------------
Nt = ceil(Tfinal / dt);
time = 0;
for n = 1:Nt
    % Explicit transport step
    h_transport = zeros(size(h));
    for j = 1:Nv
        vj = v_nodes(j);
        if vj > 0
            h_shift = [h(j,end), h(j,1:end-1)];
            dhdx = (h(j,:) - h_shift) / dx;
        elseif vj < 0
            h_shift = [h(j,2:end), h(j,1)];
            dhdx = (h_shift - h(j,:)) / dx;
        else
            dhdx = zeros(1,Nx);
        end
        h_transport(j,:) = h(j,:) - dt * vj * dhdx;
    end

    % Implicit collision step
    h_new = zeros(size(h));
    for ix = 1:Nx
        rhs = h_transport(:,ix);
        y = Plu * rhs;
        z = Llu \ y;
        h_new(:,ix) = Ulu \ z;
    end
    h = h_new;
    time = time + dt;

    % Diagnostics every 10% of steps
    if mod(n, round(Nt/10))==0 || n==1
        mid = round(Nx/2);
        f_micro_mid = Phi_plot * h(:,mid);

        % macro moments from micro
        moment_mid = W' * h(:,mid);
        a_mid = M_inv * moment_mid;
        Phi_macro_plot = [ones(length(vq_plot),1), vq_plot, vq_plot.^2];
        f_macro_mid = Phi_macro_plot * a_mid;
        f_total_mid = f_macro_mid + f_micro_mid;

        figure(1); clf;
        subplot(2,1,1);
        plot(vq_plot, f_total_mid,'k-','LineWidth',1.4); hold on;
        plot(v_nodes, h(:,mid),'ro','MarkerFaceColor','r');
        plot(vq_plot, f_macro_mid,'b--','LineWidth',1.2);
        legend('f_{total}','h (nodes)','f_{macro}');
        title(sprintf('Mid cell, time=%.4f', time)); xlabel('v'); ylabel('f(v)'); grid on;

        subplot(2,1,2);
        plot(x, W' * h,'LineWidth',1.4); xlabel('x'); ylabel('Macro moments'); title('Micro moments over x'); grid on;
        drawnow;
    end
end

fprintf('Finished. time ~ %.4f\n', time);
