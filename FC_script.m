clear;clc;
close all

%
%   A script that shows how to construct an FC extension using the FC-Gram
%   algorithm.
%
% Author:
%
%   Daniel Leibovici
%   email: dleibovi@caltech.edu
%

%% Set up the domain and FC data structures
n = 500; % number of points in the initial domain
x_a = 0; x_b = 1; % The beginning and end of the Cartesian grid
h = (x_b - x_a)/(n-1);
x = linspace(x_a, x_b, n).';

d =  5; % Number of Gram polynomial interpolation points
C = 25; % Number of continuation points
E = C;
Z = 12;


fourPts = n + C; % Number of points in the extended grid
prd = fourPts*h; % Extended period
if (mod(fourPts, 2) == 0)
    k = transpose([0:fourPts/2, -fourPts/2+1:-1]);
else
    k = transpose([0:(fourPts-1)/2, -(fourPts-1)/2:-1]);
end

% Loading continuation matrices
if(exist(['FC_data/A_d',num2str(d),'_C_', num2str(C), '.mat']) == 0 || ...
   exist(['FC_data/Q_d',num2str(d),'_C_', num2str(C), '.mat']) == 0 || ...
   exist(['FC_data/Q_tilde_d',num2str(d),'_C_', num2str(C), '.mat']) == 0)
    disp('FC data not found. Generating FC operators... \n');
    generate_bdry_continuations(d, C, E, Z, n_over, num_digits);
end
load(['FC_data/A_d',num2str(d),'_C_', num2str(C), '.mat']);
load(['FC_data/Q_d',num2str(d),'_C_', num2str(C), '.mat']);
load(['FC_data/Q_tilde_d',num2str(d),'_C_', num2str(C), '.mat']);

% Building matrices used to produce the continuation
[ArQr, AlQl, ArQ_tilder, AlQ_tildel] = build_cont_mat(A, Q, Q_tilde);


% Here, enter the desired function; for example y(x) = exp(cos(x))
y =  x.^2;


% Produce the FC extension, assuming Dirichlet boundary conditions:
BC = [0, 0; 0, 0]; % Dirichlet boundary condition, i.e. we don't impose
                   % anything on the derivative.
[~, ycont] = fcont_gram_blend(y, d, C, ArQr, AlQl, BC);
xcont = [x; x_b + h*(1 : C).'];

figure
plot(x, y, 'b', xcont, ycont, 'r--');
legend('Initial data', 'FC extended data')
