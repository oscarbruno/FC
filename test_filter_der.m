clear all
close all



d =  5; % Number of Gram polynomial interpolation points
C = 25; % Number of continuation points
E = C;
Z = 12;

load(['FC_data/A_d',num2str(d),'_C_', num2str(C), '.mat']);
load(['FC_data/Q_d',num2str(d),'_C_', num2str(C), '.mat']);
load(['FC_data/Q_tilde_d',num2str(d),'_C_', num2str(C), '.mat']);

% Pre-building matrices used to produce the continuation (for increased
% performance)
[AQ, AlQl, ArQ_tilder, AlQ_tildel] = build_cont_mat(A, Q, Q_tilde);


x_a = 0; x_b = 1; % The beginning and end of the Cartesian grid

N_test = [10; 20; 50; 100; 200; 300; 400; 600; 800; 1000; 5000; 10000];
u = @(x) exp(sin(x));
% du = @(x) 5*cos(5*x).*exp(sin(5*x)); 
du = @(x) (cos(x).^2 - sin(x)).*exp(sin(x)); 
% u = @(x) cos(x);
% du = @(x) -cos(x); 
BC = [0, 0; 0, 0];

error = zeros(size(N_test));
error_filt = zeros(size(N_test));
x_fine = linspace(x_a, x_b, N_test(end)).';

for i = 1 : length(N_test)
    n = N_test(i); % number of points in the initial domain
    h = (x_b - x_a)/(n-1);
    x = linspace(x_a, x_b, n).';


    fourPts = n + C; % Number of points in the extended grid
    prd = fourPts*h; % Extended period
    if (mod(fourPts, 2) == 0)
        k = transpose([0:fourPts/2, -fourPts/2+1:-1]);
    else
        k = transpose([0:(fourPts-1)/2, -(fourPts-1)/2:-1]);
    end

   C0 = exp(1i * 2 * pi /prd * k * x_fine.');


    alpha = 16*log(10);
    p = 3*fourPts/5;

%     alpha = 16*log(10);
%     p = 50;
%     alpha = 10;
%     p = 7;
    filter_coeffs = exp(- alpha * (2*k/fourPts).^(2*p));
    % Derivative spectral coefficients
    der_coeffs = 1i * 2*pi / prd * k;    
    der_coeffs2 = (1i * 2*pi / prd * k).^2; 
    
    [u_der, u_der_coeffs] = fc_der(u(x), der_coeffs2, 1, d, C, AQ, AlQl, BC, h);
    u_der_filt = fc_der(u(x), der_coeffs2, filter_coeffs, d, C, AQ, AlQl, BC, h); 

%     u_der = C0.' * (u_der_coeffs .* der_coeffs);
%     u_der_filter = C0.' * (u_der_coeffs .* der_coeffs .* filter_coeffs);
%     
%     u_der = real(u_der(1 : N_test(end)));
%     u_der_filter = real(u_der_filter(1 : N_test(end)));
    
    
%     error(i) = max(abs(du(x_fine) - u_der));
%     error_filt(i) = max(abs(du(x_fine) - u_der_filter));
    error(i) = max(abs(du(x) - u_der));
    error_filt(i) = max(abs(du(x) - u_der_filt));
end

% semilogy(N_test, error, 'b', N_test, error_filt, 'r');
plot(N_test, log(error), 'b', N_test, log(error_filt), 'r');


