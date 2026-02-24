clc; clear; addpath('utils');

%% Simulation settings
c = 0.1:0.2:0.9;      % c=c_s+c_f,
                      % c_s is a ratio of scattering c.s. to total c.s.
                      % c_f is a ratio of fission c.s. to total c.s.,
                      % where 'c.s.' represents 'cross section'.
r = 0;                % truncation order of Legendre expansion
n = 30;               % number of quadrature points
fs = 1;               % interaction probability
L = 5;                % sptial length
N = 200;              % number of grid points
x = linspace(0,L,N);  % grid points

%% Solve by inverse Fourier transform
phi0 = zeros(numel(c),numel(x));
for j = 1:numel(c)
    fprintf('Isotropic source test with c=%.1f\n',c(j));
    for i = 1:numel(x)
        fprintf('=> draw at x=%.4f\n',x(i));
        g = @(k) fn_phi0(k,r,n,c(j),fs);
        phi0(j,i) = cifft(g,x(i),10,0.45);
    end
end

%% Visulization
figure(1)
semilogy(x(2:end),phi0(1,2:end),'-', ...
    'LineWidth',1,'DisplayName',sprintf('c=%.1f',c(1))); 
hold on;
for j = 2:numel(c)
    semilogy(x(2:end),phi0(j,2:end),'--', ...
        'LineWidth',1,'DisplayName',sprintf('c=%.1f',c(j)));
end
hold off;
title('Isotropic source')
xlabel('$x$','Interpreter','latex','FontSize',14);
ylabel('Scalar Flux','FontSize',14);
xlim([-1 L+1]);
ylim([1e-4 1e1]);
legend('Location','best');
save('data_og1d_plane_isotropic_source_vc.mat');

%% Utils
function phi0 = fn_phi0(k,r,n,c,fs)
    % compute matrix T and A
    T = zeros(r+1,r+1);
    for m = 0:r
        for l = 0:r
            T(m+1,l+1) = gaussint( ...
                @(mu)legdval(m,mu).*legdval(l,mu)./(1+1j*k*mu),-1,1,n);
        end
    end
    T = T/2;
    A = real(1j.^((0:r)-(0:r)').*T); % A is a real-valued matrix
    % compute matrix Q and S
    Q = zeros(1,r+1);
    for m = 0:r
        Q(1,m+1) = gaussint( ...
                @(mu)legdval(m,mu)./(1+1j*k*mu)/2,-1,1,n);
    end
    S = 1;
    % compute phi0
    R = 1j.^(-(0:r)).*Q;
    L = eye(r+1,r+1)-(2*(0:r)+1).*c.*fs.*A;
    X = (L\real(R(:)))+1j*(L\imag(R(:)));
    phi = (1j.^(0:r))'.*X*S;
    phi0 = phi(1);
end


