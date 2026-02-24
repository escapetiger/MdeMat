%% Function: continuous inverse fast Fourier transform
%  Args:
%   g - function in Fourier space;
%   x - physical coordinates;
%   n - truncation order;
%   pe - parameter for Euler-Knopp acceleration;
function f = cifft(g,x,n,pe)
    if nargin<4
        pe = 0.45;
    end
    a = @(j)integral_segment(g,j,x);
    if pe == 0 % without EK acceleration
        f = 0;
        for k = 0:n
           f = f+(-1).^(k).*a(k);
        end
    else % with EK accleration
        A = (1+pe).^(1:n+1);
        B = binom(0:n,0:n);
        U = pe.^(0:n);
        V = zeros(1,n+1);
        for k = 0:n
            V(k+1) = (-1)^k*a(k);
        end
        f = 0;
        for k = 0:n
            c = B(k+1,1:k+1)./A(k+1);
            u = U(1:k+1);
            v = V(k+1:-1:1);
            f = f+sum(c.*u.*v);
        end
    end
end

function H = integrand(g,k,u)
    G = g(k); GR = real(G); GI = imag(G);
    H = GR.*cos(u)-GI.*sin(u);
end

function I = integral_segment(g,j,x,n_gauss)
    if nargin<4
        n_gauss = 40;
    end
    if x==0
        T = @(t) (1+t)./(1-t);
        I = gaussint(@(t)real(g(T(t)))./(1-t).^2,-1,1,n_gauss)*2/pi;
    elseif x>0
        if x<=0.1 && j == 0
            I = gaussint(@(k)integrand(g,k,k*x),0,pi/x,n_gauss)/pi;
        else 
            T = @(u,j)(u+j*pi)./x;
            I = gaussint(@(u)integrand(g,T(u,j),u),0,pi,n_gauss)./(pi*x);
        end
    elseif x<0 
        if x>=-0.1 && j==0
            I = -gaussint(@(k)integrand(g,k,k*x),-pi/x,0,n_gauss)/pi;
        else
            T = @(u,j)(u-j*pi)./x;
            I = -gaussint(@(u)integrand(g,T(u,j),u),-pi,0,n_gauss)./(pi*x);
        end
    end
end
