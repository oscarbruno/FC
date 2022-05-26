clear all
close all


%% Set up the domain and FC data structures

n = 100; % number of points in the initial domain
x_a = 0; x_b = 1; % The beginning and end of the Cartesian grid
h = (x_b - x_a)/(n-1);
x = linspace(x_a, x_b, n).';

d =  5; % Number of Gram polynomial interpolation points
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
if(exist(['FC_data/A_d',num2str(d),'_C_', num2str(C), '.mat']) == 0 | ...
   exist(['FC_data/Q_d',num2str(d),'_C_', num2str(C), '.mat']) == 0 | ...
   exist(['FC_data/Q_tilde_d',num2str(d),'_C_', num2str(C), '.mat']) == 0)
    disp('FC data not found. Generating FC operators... \n');
    generate_bdry_continuations(d, C, E, Z, n_over, num_digits);
end
load(['FC_data/A_d',num2str(d),'_C_', num2str(C), '.mat']);
load(['FC_data/Q_d',num2str(d),'_C_', num2str(C), '.mat']);
load(['FC_data/Q_tilde_d',num2str(d),'_C_', num2str(C), '.mat']);

% Pre-building matrices used to produce the continuation (for increased
% performance)
[ArQr, AlQl, ArQ_tilder, AlQ_tildel] = build_cont_mat(A, Q, Q_tilde);

% Derivative spectral coefficients
der_coeffs = 1i* 2*pi / prd * k;
der_coeffs2 = der_coeffs.^2;

% Filter spectral coefficients
alpha = 35;
p = 10;
filter_coeffs = 1;

%% A 2nd order Wave Equation example

deltat = h/32; % Time step for the AB4 scheme
maxit = ceil(2/deltat); % Number of iterations
T = maxit*deltat; % Final computation time

u = zeros(length(x), maxit); % Storing the solution
g = zeros(length(x), maxit); % Storing the solution's time derivative


uexact = @(x, t) cos(t).*cos(x); % The exact solution
dtuexact = @(x, t) -sin(t).*cos(x); % The time derivative of uexact
dxuexact = @(x, t) - cos(t).*sin(x); % The spatial derivative of uexact


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
    
%     BC0 = [1, 1; dxuexact(x_a, t), dxuexact(x_b, t)];
%     BC1 = [1, 1; dxuexact(x_a, t - deltat), dxuexact(x_b, t - deltat)];
%     BC2 = [1, 1; dxuexact(x_a, t - 2*deltat), dxuexact(x_b, t - 2*deltat)];
%     BC3 = [1, 1; dxuexact(x_a, t - 3*deltat), dxuexact(x_b, t - 3*deltat)];
% 
%     g(:, it + 1) = g(:, it) + ... 
%         deltat/24*(55*fc_der(u(:, it), der_coeffs2, filter_coeffs, d, C, ArQ_tilder, AlQ_tildel, BC0, h) ...
%         - 59*fc_der(u(:, it - 1), der_coeffs2, filter_coeffs, d, C, ArQ_tilder, AlQ_tildel, BC1, h)...
%         + 37*fc_der(u(:, it - 2), der_coeffs2, filter_coeffs, d, C, ArQ_tilder, AlQ_tildel, BC2, h)...
%         - 9*fc_der(u(:, it - 3), der_coeffs2, filter_coeffs, d, C, ArQ_tilder, AlQ_tildel, BC3, h)); 
%     
%     u(:, it + 1) = u(:, it) + deltat/24*(55*g(:, it) - 59*g(:, it - 1) + ...
%                                         37*g(:, it - 2) - 9*g(:, it - 3)); 
%     BC = [1, 1; dxuexact(x_a, t + deltat), dxuexact(x_b, t + deltat)];                                
%     u(:, it + 1) = fc_der(u(:, it + 1), 1, filter_coeffs, d, C, ArQ_tilder, AlQ_tildel, BC, h);    


    % Left Dirichlet and Right Neumann boundary conditions
    
    BC0 = [0, 1; dxuexact(0, t), dxuexact(1, t)];
    BC1 = [0, 1; dxuexact(0, t - deltat), dxuexact(1, t - deltat)];
    BC2 = [0, 1; dxuexact(0, t - 2*deltat), dxuexact(1, t - 2*deltat)];
    BC3 = [0, 1; dxuexact(0, t - 3*deltat), dxuexact(1, t - 3*deltat)];
    g(:, it + 1) = g(:, it) + deltat/24*(...
          55*fc_der(u(:, it), der_coeffs2, filter_coeffs, d, C, ArQ_tilder, AlQl, BC0, h) ...
        - 59*fc_der(u(:, it - 1), der_coeffs2, filter_coeffs, d, C, ArQ_tilder, AlQl, BC1, h)...
        + 37*fc_der(u(:, it - 2), der_coeffs2, filter_coeffs, d, C, ArQ_tilder, AlQl, BC2, h)...
        - 9*fc_der(u(:, it - 3), der_coeffs2, filter_coeffs, d, C, ArQ_tilder, AlQl, BC3, h));
    u(:, it + 1) = u(:, it) + deltat/24*(55*g(:, it) - 59*g(:, it - 1) + ...
                                        37*g(:, it - 2) - 9*g(:, it - 3));
    BC = [0, 1; dxuexact(0, t), dxuexact(1, t + deltat)];
    u(:, it + 1) = fc_der(u(:, it + 1), 1, filter_coeffs, d, C, ArQ_tilder, AlQl, BC, h);
    u(1, it + 1) = uexact(0, t + deltat); 
    
    

%     Left and Right Dirichlet boundary conditions

%     BC0 = [0, 0; dxuexact(0, t), dxuexact(1, t)];
%     BC1 = [0, 0; dxuexact(0, t - deltat), dxuexact(1, t - deltat)];
%     BC2 = [0, 0; dxuexact(0, t - 2*deltat), dxuexact(1, t - 2*deltat)];
%     BC3 = [0, 0; dxuexact(0, t - 3*deltat), dxuexact(1, t - 3*deltat)];
% 
%     g(:, it + 1) = g(:, it) + deltat/24*(...
%           55*fc_der(u(:, it), 2, filter, prd, k, d, C, A, Q, Q, BC0) ...
%         - 59*fc_der(u(:, it - 1), 2, filter, prd, k, d, C, A, Q, Q, BC1)...
%         + 37*fc_der(u(:, it - 2), 2, filter, prd, k, d, C, A, Q, Q, BC2)...
%         - 9*fc_der(u(:, it - 3), 2, filter, prd, k, d, C, A, Q, Q, BC3));    
%     u(:, it + 1) = u(:, it) + deltat/24*(55*g(:, it) - 59*g(:, it - 1) + ...
%                                         37*g(:, it - 2) - 9*g(:, it - 3));
%     u(1, it + 1) = uexact(0, t + deltat);
%     u(end, it + 1) = uexact(1, t + deltat);
    
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
err = u(:, it) - uexact(x, plottime);
semilogy(x, abs(err));
title(['Error at t = ', num2str(plottime)]);
fprintf('Maximum error: %1.3e\n', max(abs(err)));
