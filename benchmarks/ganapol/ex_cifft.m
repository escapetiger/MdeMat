clc; clear; addpath('utils');

%% alternating series
% N = 50;
% c = 0.99;
% a = @(k) c.^k;
% S = sum((-1).^(0:N).*a(0:N));
% Se = euler_knopp(a,N,0.5);
% T = 1/(1+c);

%% Fourier inversion: Test 1
n1 = 10; n2 = 30000;
x = 0.1;
g1 = @(k)-1j*pi*k./(1+k.^2);
h1 = @(k)k./(1+k.^2).*sin(k*x);
I1 = cifft(g1,x,n1,0.45);
J1 = integral(h1,0,n2*pi);

%% Fourier inversion: Test 2
n1 = 10; n2 = 100;
x = 1;
g2 = @(k)1j*pi*sign(k)./(abs(k)+pi).^2;
h2 = @(k)sin(k*x)./k.^2;
I2 = cifft(g2,x,n1,0.45);
J2 = integral(h2,pi,n2*pi);

%% Fourier inversion: Test 3
n1 = 10; n2 = 100;
x = 1;
g3 = @(k)-1j*pi*exp(-abs(k)./2)./k;
h3 = @(k)sin(k*x)./(k).*exp(-k/2);
I3 = cifft(g3,x,n1,0.45);
J3 = integral(h3,0,n2*pi);

%% Fourier inversion: Test 4
n1 = 10; n2 = 100;
x = 4;
g4 = @(k)2*pi*exp(-2*abs(k))./(1+exp(-4*abs(k)));
h4 = @(k)sech(2*k).*cos(k*x);
I4 = cifft(g4,x,n1,0.45);
J4 = integral(h4,0,n2*pi);

%% Fourier inversion: Test 5
% NOTE: J is not converged
n1 = 10; n2 = 1000;
x = 1;
g5 = @(k)-1j*pi*sign(k)./(abs(k).^0.1);
h5 = @(k)sin(k*x)./(k.^0.1);
I5 = cifft(g5,x,n1,0.45);
J5 = integral(h5,0,n2*pi);

%% Fourier inversion: Test 6
% NOTE: J is not converged
n1 = 10; n2 = 1000;
x = 1;
g6 = @(k)-1j*pi./(2*k.^(3/2)).*sqrt(abs(k));
h6 = @(k)sqrt(k).*sin(100*x*k)/2;
I6 = cifft(g6,x,n1,0.45);
J6 = integral(h6,0,n2*pi);

%% Fourier inversion: Test 7
n1 = 10; n2 = 100;
x = 1;
g7 = @(k)pi./(k.^2+4);
h7 = @(k)cos(x*k)./(k.^2+4);
I7 = cifft(g7,x,n1,0.45);
J7 = integral(h7,0,n2*pi);

%% Fourier inversion: Test 8
n1 = 10; n2 = 100;
x = 1;
g8 = @(k)pi*((2+k)./(2+2*k+k.^2)+(2-k)./(2-2*k+k.^2));
h8 = @(k)(2+k)./(2+2*k+k.^2).*cos(k*x);
I8 = cifft(g8,x,n1,0.45);
J8 = integral(h8,-n2*pi,n2*pi);

%% Fourier inversion: Test 9
n1 = 5; n2 = 10;
x = 10;
g9 = @(k)-1j*pi*exp(-(abs(k)+1/abs(k))).*sign(k)./sqrt(abs(k));
h9 = @(k)exp(-(k+1./k))./sqrt(k).*sin(k*x);
I9 = cifft(g9,x,n1,0.45);
J9 = integral(h9,0,n2*pi);