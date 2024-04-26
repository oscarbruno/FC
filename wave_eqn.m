clear all
close all

%
%   A 1D Second order wave equation script.
%
% Author:
%
%   Daniel Leibovici
%   email: dleibovi@caltech.edu
%

%% Set up the domain and FC data structures

n = 801; % number of points in the initial domain
x_a = 0; x_b = 1; % The beginning and end of the Cartesian grid
h = (x_b - x_a)/(n-1);
x = linspace(x_a, x_b, n).';

d =  5; % Number of Gram polynomial interpolation points
C = 27; % Number of continuation points
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
if(exist(['FC_data/A_d',num2str(d),'_C', num2str(C), '.mat']) == 0 | ...
   exist(['FC_data/Q_d',num2str(d),'_C', num2str(C), '.mat']) == 0 | ...
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
der_coeffs2 = der_coeffs.^2;

% Filter spectral coefficients
alpha = 35;
p = 10;
% filter_coeffs = 1;
% alpha = 16 * log(10);
% p = 3 * fourPts / 5;
% filter_coeffs = exp(- alpha * (2 * k / fourPts) .^ (2 * p));


%% A 2nd order Wave Equation example

deltat = h / 32; % Time step for the AB4 scheme


cmax = 1;
p = 4;
alpha = - cmax * deltat / h * log(0.01);
filter_coeffs = exp(- alpha * (2 * k / fourPts) .^ (2 * p));

T = 2; % Final computation time
maxit = T / deltat; % Number of iterations

u = zeros(length(x), maxit); % Storing the solution
g = zeros(length(x), maxit); % Storing the solution's time derivative


% uexact = @(x, t) cos(t).*cos(x); % The exact solution
% dtuexact = @(x, t) -sin(t).*cos(x); % The time derivative of uexact
% dxuexact = @(x, t) - cos(t).*sin(x); % The spatial derivative of uexact
uexact = @(x, t) exp(-300 * (x - t + 0.5).^2); % The exact solution
dtuexact = @(x, t) 600 * (x - t + 0.5) .* exp(-300 * (x - t + 0.5).^2); % The time derivative of uexact
dxuexact = @(x, t) - 600 * (x - t + 0.5) .* exp(-300 * (x - t + 0.5).^2); % The spatial derivative of uexact


u0x = uexact(x, 0); % The initial condition

% Pre-setting the initial condition and the first time steps to equal the
% exact solution
u(:, 1) = u0x; 
u(:, 2) = uexact(x, deltat);
u(:, 3) = uexact(x, 2*deltat); 
u(:, 4) = uexact(x, 3*deltat); 

% Pre-setting the initial condition of g and the first time steps to equal
% the time derivative of the exact solution
g(:, 1) = dtuexact(x, 0);
g(:, 2) = dtuexact(x, deltat);
g(:, 3) = dtuexact(x, 2*deltat);
g(:, 4) = dtuexact(x, 3*deltat);

tic;

t = 3*deltat;
it = 4;

% Time stepping with the AB4 scheme
while (it <= maxit)

    t
    
    % Left and right Neumann boundary conditions
    

    BC0 = [1, 1; dxuexact(x_a, t), dxuexact(x_b, t)];
    BC1 = [1, 1; dxuexact(x_a, t - deltat), dxuexact(x_b, t - deltat)];
    BC2 = [1, 1; dxuexact(x_a, t - 2*deltat), dxuexact(x_b, t - 2*deltat)];
    BC3 = [1, 1; dxuexact(x_a, t - 3*deltat), dxuexact(x_b, t - 3*deltat)];

    g(:, it + 1) = g(:, it) + ... 
        deltat/24*(55*fc_der(u(:, it), der_coeffs2, filter_coeffs, d, ...
              C, A, Q_tilde, Q_tilde, BC0, h) ...
        - 59*fc_der(u(:, it - 1), der_coeffs2, filter_coeffs, d, C, ... 
              A, Q_tilde, Q_tilde, BC1, h)...
        + 37*fc_der(u(:, it - 2), der_coeffs2, filter_coeffs, d, C, ... 
              A, Q_tilde, Q_tilde, BC2, h)...
        - 9*fc_der(u(:, it - 3), der_coeffs2, filter_coeffs, d, C, ... 
              A, Q_tilde, Q_tilde, BC3, h)); 
    
    u(:, it + 1) = u(:, it) + deltat/24*(55*g(:, it) ... 
           - 59*g(:, it - 1) + 37*g(:, it - 2) - 9*g(:, it - 3)); 
    BC = [1, 1; dxuexact(x_a, t + deltat), dxuexact(x_b, t + deltat)];                                
    u(:, it + 1) = fc_der(u(:, it + 1), 1, filter_coeffs, d, C, A, ... 
          Q_tilde, Q_tilde, BC, h);   
      



    % Left Neumann and right Dirichlet boundary conditions
    

%     BC0 = [1, 0; dxuexact(x_a, t), dxuexact(x_b, t)];
%     BC1 = [1, 0; dxuexact(x_a, t - deltat), dxuexact(x_b, t - deltat)];
%     BC2 = [1, 0; dxuexact(x_a, t - 2*deltat), dxuexact(x_b, t - 2*deltat)];
%     BC3 = [1, 0; dxuexact(x_a, t - 3*deltat), dxuexact(x_b, t - 3*deltat)];
% 
%     g(:, it + 1) = g(:, it) + ... 
%         deltat/24*(55*fc_der(u(:, it), der_coeffs2, filter_coeffs, d, ...
%               C, A, Q_tilde, Q, BC0, h) ...
%         - 59*fc_der(u(:, it - 1), der_coeffs2, filter_coeffs, d, C, ... 
%               A, Q_tilde, Q, BC1, h)...
%         + 37*fc_der(u(:, it - 2), der_coeffs2, filter_coeffs, d, C, ... 
%               A, Q_tilde, Q, BC2, h)...
%         - 9*fc_der(u(:, it - 3), der_coeffs2, filter_coeffs, d, C, ... 
%               A, Q_tilde, Q, BC3, h)); 
% 
%     
%     u(:, it + 1) = u(:, it) + deltat/24*(55*g(:, it) ... 
%            - 59*g(:, it - 1) + 37*g(:, it - 2) - 9*g(:, it - 3)); 
%     BC = [0, 1; dxuexact(x_a, t + deltat), dxuexact(x_b, t + deltat)];                                
%     u(:, it + 1) = fc_der(u(:, it + 1), 1, filter_coeffs, d, C, A, ... 
%           Q_tilde, Q, BC, h);  
%     u(end, it + 1) = uexact(1, t + deltat); 


    % Left Dirichlet and Right Neumann boundary conditions
    
%     BC0 = [0, 1; dxuexact(0, t), dxuexact(1, t)];
%     BC1 = [0, 1; dxuexact(0, t - deltat), dxuexact(1, t - deltat)];
%     BC2 = [0, 1; dxuexact(0, t - 2*deltat), dxuexact(1, t - 2*deltat)];
%     BC3 = [0, 1; dxuexact(0, t - 3*deltat), dxuexact(1, t - 3*deltat)];
%     g(:, it + 1) = g(:, it) + deltat/24*(...
%           55*fc_der(u(:, it), der_coeffs2, filter_coeffs, d, C, A, ... 
%               Q, Q_tilde, BC0, h) ...
%         - 59*fc_der(u(:, it - 1), der_coeffs2, filter_coeffs, d, C, ... 
%               A, Q, Q_tilde, BC1, h)...
%         + 37*fc_der(u(:, it - 2), der_coeffs2, filter_coeffs, d, C, ... 
%               A, Q, Q_tilde, BC2, h)...
%         - 9*fc_der(u(:, it - 3), der_coeffs2, filter_coeffs, d, C, A, ...
%               Q, Q_tilde, BC3, h));
%     u(:, it + 1) = u(:, it) + deltat/24*(55*g(:, it) ...
%           - 59*g(:, it - 1) + 37*g(:, it - 2) - 9*g(:, it - 3));
%     BC = [0, 1; dxuexact(x_a, t + deltat), dxuexact(x_b, t + deltat)];
%     u(:, it + 1) = fc_der(u(:, it + 1), 1, filter_coeffs, d, C, A, ... 
%           Q, Q_tilde, BC, h);
%     u(1, it + 1) = uexact(0, t + deltat); 
    
    

%     Left and Right Dirichlet boundary conditions

%     BC0 = [0, 0; dxuexact(0, t), dxuexact(1, t)];
%     BC1 = [0, 0; dxuexact(0, t - deltat), dxuexact(1, t - deltat)];
%     BC2 = [0, 0; dxuexact(0, t - 2*deltat), dxuexact(1, t - 2*deltat)];
%     BC3 = [0, 0; dxuexact(0, t - 3*deltat), dxuexact(1, t - 3*deltat)];
% 
%     g(:, it + 1) = g(:, it) + deltat/24*(...
%           55*fc_der(u(:, it), der_coeffs2, filter_coeffs, d, C, A, Q, ... 
%                 Q, BC0, h) ...
%         - 59*fc_der(u(:, it - 1), der_coeffs2, filter_coeffs, d, C, A, ... 
%                 Q, Q, BC1, h)...
%         + 37*fc_der(u(:, it - 2), der_coeffs2, filter_coeffs, d, C, A, ... 
%                 Q, Q, BC2, h)...
%         - 9*fc_der(u(:, it - 3), der_coeffs2, filter_coeffs, d, C, A, ... 
%                 Q, Q, BC3, h));
%     u(:, it + 1) = u(:, it) + deltat/24*(55*g(:, it) - ... 
%             59*g(:, it - 1) + 37*g(:, it - 2) - 9*g(:, it - 3));
%     u(1, it + 1) = uexact(0, t + deltat);
%     u(end, it + 1) = uexact(1, t + deltat);
%     
    t = t + deltat;
    it = it + 1;

end


toc;

plottime = T;

figure
plot(x, u(:, it));
title(['Solution at t = ', num2str(plottime)]);


figure(1)
plot(x, uexact(x, T));
title(['exact solution at t = ', num2str(plottime)]);

figure
err = u(:, it) - uexact(x, t);
semilogy(x, abs(err));
title(['Error at t = ', num2str(plottime)]);
fprintf('Maximum error: %1.3e\n', max(abs(err)));
