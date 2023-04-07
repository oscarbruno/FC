clear all
close all

%
%   A test for the accuracy of the FC approximation.
%
% Author:
%
%   Daniel Leibovici
%   email: dleibovi@caltech.edu
%

%% Set up the domain and FC data structures

n = 1000; % number of points in the initial domain
x_a = 0; x_b = 1; % The beginning and end of the Cartesian grid
h = (x_b - x_a)/(n-1);
x = linspace(x_a, x_b, n).';

d =  10; % Number of Gram polynomial interpolation points
C = 27; % Number of continuation points

Z = 12;
E = C;
n_over = 20;
num_digits = 256;

fourPts = n + C; % Number of points in the extended grid
prd = fourPts*h; % Extended period
if (mod(fourPts, 2) == 0)
    k = transpose([0:fourPts/2, -fourPts/2+1:-1]);
else
    k = transpose([0:(fourPts-1)/2, -(fourPts-1)/2:-1]);
end


% Loading continuation matrices
if(exist(['FC_data/A_d',num2str(d),'_C', num2str(C), '.mat']) == 0 | ...
   exist(['FC_data/Q_d',num2str(d),'_C', num2str(C), '.mat']) == 0 | ...
   exist(['FC_data/Q_tilde_d',num2str(d),'_C', num2str(C), '.mat']) == 0)
    disp('FC data not found. Generating FC operators... \n');
    generate_bdry_continuations(d, C, E, Z, n_over, num_digits);
end

% generate_bdry_continuations(d, C, E, Z, n_over, num_digits);
load(['FC_data/A_d',num2str(d),'_C', num2str(C), '.mat']);
load(['FC_data/Q_d',num2str(d),'_C', num2str(C), '.mat']);
load(['FC_data/Q_tilde_d',num2str(d),'_C', num2str(C), '.mat']);


% A = dlmread(['FC_data/A', num2str(d), 'C', num2str(C), '.dat']);
% Q = dlmread(['FC_data/Q', num2str(d), 'C', num2str(C), '.dat']);
% Q_tilde = dlmread(['FC_data/Q_tilde', num2str(d), 'C', num2str(C), '.dat']);
% AlQl = dlmread(['FC_data/AlQl', num2str(d), 'C', num2str(C), '.dat']);
% ArQr = dlmread(['FC_data/ArQr', num2str(d), 'C', num2str(C), '.dat']);
% AlQ_tildel = dlmread(['FC_data/AlQ_tildel', num2str(d), 'C', num2str(C), '.dat']);
% ArQ_tilder = dlmread(['FC_data/ArQ_tilder', num2str(d), 'C', num2str(C), '.dat']);

% A = vpa(dlmread(['FC_data/A', num2str(d), 'C', num2str(C), '.dat']), num_digits);
% Q = vpa(dlmread(['FC_data/Q', num2str(d), 'C', num2str(C), '.dat']), num_digits);
% Q_tilde = vpa(dlmread(['FC_data/Q_tilde', num2str(d), 'C', num2str(C), '.dat']), num_digits);
% AlQl = vpa(dlmread(['FC_data/AlQl', num2str(d), 'C', num2str(C), '.dat']), num_digits);
% ArQr = vpa(dlmread(['FC_data/ArQr', num2str(d), 'C', num2str(C), '.dat']), num_digits);
% AlQ_tildel = vpa(dlmread(['FC_data/AlQ_tildel', num2str(d), 'C', num2str(C), '.dat']), num_digits);
% ArQ_tilder = vpa(dlmread(['FC_data/ArQ_tilder', num2str(d), 'C', num2str(C), '.dat']), num_digits);


% Pre-building matrices used to produce the continuation (for increased
% performance)
[ArQr, AlQl, ArQ_tilder, AlQ_tildel] = build_cont_mat(A, Q, Q_tilde);
A = double(A);
Q = double(Q);
Q_tilde = double(Q_tilde);
ArQr = double(ArQr);
AlQl = double(AlQl);
ArQ_tilder = double(ArQ_tilder);
AlQ_tildel = double(AlQ_tildel);


% Derivative spectral coefficients
der_coeffs0 = 1;
der_coeffs = 1i* 2 * pi / prd * k;
der_coeffs2 = der_coeffs.^2;

% Filter spectral coefficients
alpha = 35;
p = 10;
filter_coeffs = 1;



%% Evaluate the function and perform continuations
%f = @(y) y.^2 - 1*y + 1;
%fder = @(y) 2*y - 1;
%f2der = @(y) 2;

a = 5.4*pi;
b = 2.7*pi;
c = 2*pi;
f = @(y) exp(sin(a*y - b) - cos(c*y));
% f = @(y) exp(cos(2 * pi * y));
fder = @(y) (a*cos(a*y - b) + c*sin(c*y)).*exp(sin(a*y - b) - cos(c*y));
f2der = @(y) -(a^2*sin(a*y - b) - c^2*cos(c*y) - (a*cos(a*y - b) + c*sin(c*y)).^2) .* exp(sin(a*y - b) - cos(c*y));

fx = f(x);
BC0 = [0, 0; 0, 0];

[~, fc_coeffs] = fc_der(fx, der_coeffs0, filter_coeffs, d, C, A, Q, Q, ArQr, AlQl, BC0, h);

z = linspace(x_a, x_b, 20000);

fc_fine_grid = real(exp(2*pi*1i*(z - x_a).'*k.'/prd) * ...
        fc_coeffs(:));
ferr = abs(fc_fine_grid - f(z).');
fprintf('||f(x)||: %1.3e\n', max(ferr));

figure
plot(z, log10(ferr));
xlabel('x', 'Fontsize', 15);
ylabel('log_1_0(e)', 'Fontsize', 15);







