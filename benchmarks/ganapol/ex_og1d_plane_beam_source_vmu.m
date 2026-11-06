function result = ex_og1d_plane_beam_source_vmu(options)
arguments
    options.Plot (1, 1) logical = true
    options.OutputFile (1, 1) string = ""
    options.MaxPoints (1, 1) double {mustBePositive} = Inf
end
oldPath = path;
pathGuard = onCleanup(@() path(oldPath));
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));

%% Simulation settings
c = 0.7;              % c=c_s+c_f,
                      % c_s is a ratio of scattering c.s. to total c.s.
                      % c_f is a ratio of fission c.s. to total c.s.,
                      % where 'c.s.' represents 'cross section'.
r = 0;                % truncation order of Legendre expansion
n = 30;               % number of quadrature points
fs = 1;               % interaction probability
L = 2;                % sptial length
N = 200;              % number of grid points
x = linspace(-L,L,N); % grid points
mu0 = 0.1:0.2:0.9;    % beam direction

x = x(1:min(numel(x), options.MaxPoints));

%% Solve by inverse Fourier transform
phi0 = zeros(numel(mu0),numel(x));
for j = 1:numel(mu0)
    fprintf('Beam source test with mu0=%.1f\n',mu0(j));
    for i = 1:numel(x)
        fprintf('=> draw at x=%.4f\n',x(i));
        g = @(k) fn_phi0(k,r,n,c,fs,mu0(j));
        phi0(j,i) = cifft(g,x(i),10,0.45);
    end
end

if options.Plot
%% Visulization
figure(1)
semilogy(x,phi0(1,:),'-', ...
    'LineWidth',1,'DisplayName',sprintf('%s_0=%.1f',char(181),mu0(1))); 
hold on;
for j = 2:numel(mu0)
    semilogy(x,phi0(j,:),'--', ...
        'LineWidth',1,'DisplayName',sprintf('%s_0=%.1f',char(181),mu0(j)));
end
hold off;
title('Beam source','FontSize',16);
xlabel('$x$','Interpreter','latex','FontSize',14);
ylabel('Scalar Flux','FontSize',14);
xlim([-L-1 L+1]);
ylim([7e-2 1.2e1]);
legend('Location','best');
end

result = struct('Points', x, 'Flux', phi0, 'ScatteringRatio', c, ...
    'Phase', fs, 'BeamDirection', mu0, 'AngularPoints', n, 'FourierOrder', 10);
if options.OutputFile ~= ""
    mdemat.tools.io.writeRecord(options.OutputFile, result);
end

end


%% Utils
function phi0 = fn_phi0(k,r,n,c,fs,mu0)
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
        Q(1,m+1) = legdval(m,mu0)./(1+1j*k*mu0);
    end
    S = 1;
    % compute phi0
    R = 1j.^(-(0:r)).*Q;
    L = eye(r+1,r+1)-(2*(0:r)+1).*c.*fs.*A;
    X = (L\real(R(:)))+1j*(L\imag(R(:)));
    phi = (1j.^(0:r))'.*X*S;
    phi0 = phi(1);
end


