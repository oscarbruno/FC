clear;clc;
close all

%
%   A 1D Linear Advection equation script.
%
% Author:
%
%   Daniel Leibovici
%   email: dleibovi@caltech.edu
%

%% Set up the domain and FC data structures
n = 200; % number of points in the initial domain
x_a = 0; x_b = 1; % The beginning and end of the Cartesian grid
h = (x_b - x_a)/(n-1);
x = linspace(x_a, x_b, n).';

d =  5; % Number of Gram polynomial interpolation points
C = 25; % Number of continuation points
Z = 12;
E = C;
n_over = 20;
modes_to_reduce = 4;
num_digits = 256;

fourPts = n + C; % Number of points in the extended grid
prd = fourPts*h; % Extended period
if (mod(fourPts, 2) == 0)
    k = transpose([0:fourPts/2, -fourPts/2+1:-1]);
else
    k = transpose([0:(fourPts-1)/2, -(fourPts-1)/2:-1]);
end

% Loading continuation matrices
if(exist(['FC_data/A_d',num2str(d),'_C', num2str(C), '.mat']) == 0 || ...
   exist(['FC_data/Q_d',num2str(d),'_C', num2str(C), '.mat']) == 0 || ...
   exist(['FC_data/Q_tilde_d',num2str(d),'_C', num2str(C), '.mat']) == 0)
    disp('FC data not found. Generating FC operators... \n');
    generate_bdry_continuations(d, C, E, Z, n_over, modes_to_reduce, ...
        num_digits);
end
load(['FC_data/A_d',num2str(d),'_C', num2str(C), '.mat']);
load(['FC_data/Q_d',num2str(d),'_C', num2str(C), '.mat']);
load(['FC_data/Q_tilde_d',num2str(d),'_C', num2str(C), '.mat']);

A = double(A);
Q = double(Q);
Q_tilde = double(Q_tilde);


% Derivative spectral coefficients
der_coeffs = 1i* 2*pi / prd * k;

% Filter spectral coefficients
alpha = 10;
p = 14;
filter_coeffs = exp(- alpha * (2*k/fourPts).^p);


%% A Linear Advection Equation example

uexact = @(x, t) cos(5*pi*(x - t)); % Exact solution
u0x = uexact(x, 0); % Initial condition


tic;
deltat = 10^-6; % initially 1e-6
maxit = 1000; % initally 5000

u = u0x;
fprintf('Performing %d timesteps to t = %1.2e\n', maxit, maxit*deltat);
for it=2:maxit
    BC = [0, 0; 0, 0];
    t = (it - 1) * deltat;
    k1 = -deltat * fc_der(u, der_coeffs, filter_coeffs, d, C, A, Q, Q, ...
        BC, h);
    k2 = -deltat * fc_der(u + 1/2 * k1, der_coeffs, filter_coeffs, d, ... 
        C, A, Q, Q, BC, h);
    k3 = -deltat * fc_der(u + 1/2 * k2, der_coeffs, filter_coeffs, d, ...
        C, A, Q, Q, BC, h);
    k4 = -deltat * fc_der(u + k3, der_coeffs, filter_coeffs, d, C, A, ... 
        Q, Q, BC, h);
    u = u + 1/6*k1 + 1/3*k2 + 1/3*k3 + 1/6*k4;
    u(1) = uexact(0, t);


end
toc;

plottime = (maxit - 1)*deltat;


figure
err = u - uexact(x, plottime);
semilogy(x, abs(err));
set(gca, 'YLim', [10^-20 1]);
title(['Error at t = ', num2str(plottime)]);
fprintf('Maximum error: %1.3e\n', max(abs(err)));
