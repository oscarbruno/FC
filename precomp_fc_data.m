% This routine precomputes the Fourier Continuation matrices A_r, Q_r and
% Q_tilde_r as defined in [1]
%
% Inputs:
%   d : number of Gram interpolation points
%   C : number of continuation points
%   Z : number of zero blending points
%   E : number of extra points
%   n_over : oversampling factor
%   num_digits : number of digits used for the high precision arithmetic 
%   computation of the QR decomposition
%
% Outputs:
%   A, Q and Q_tilde, are respectively the matrices denoted A_r, Q_r in 
%   equation (26) in [1], and Q_tilde_r in equation (27) in [1]. 
%   
% Ref : 
% [1] Amlani, F., & Bruno, O. P. (2016). An FC-based spectral solver for
%     elastodynamic problems in general three-dimensional domains. 
%     Journal of Computational Physics, 307, 333-354.
%
% Author: Daniel Leibovici
% Email: dleibovi@caltech.edu
%

function [Q, Q_tilde, A] =  precomp_fc_data(d, C, Z, E, n_over, num_digits)


tol = 1e-16;
digits(num_digits);
N_coarse = d + C + Z + E;
svd_order = N_coarse;

% Forming the grids we will use
N0 = 101;
h0 = 1/(N0 - 1);
x0 = 1 - (d - 1)*h0;
h_over = h0 / n_over;
midpoint = x0 + (d + C)*h0;


interp_coarse = sym(x0 + h0 * (0 : d - 1).'); 
interp_fine = sym(x0 + h_over * linspace(0, d - 1, n_over*(d - 1) + 1).'); 
zero_fine = sym(midpoint + h_over*linspace(d + C, d - 1 + C + Z, n_over*(Z - 1) + 1).'); 
continuation = sym(x0 + h0 * (d : d - 1 + C).');


% Forming Q by using the modified Gram-Schmidt algorithm
P_coarse = sym(zeros(d, d));
P_coarse_der = sym(zeros(d, d));
for i = 0 : d - 1
    P_coarse(:, i + 1) = interp_coarse.^i;
    P_coarse_der(1 : end - 1, i + 1) = interp_coarse(1 : end - 1).^i;
    if (i == 0)
        P_coarse_der(end, i + 1) = 0;
    else
        P_coarse_der(end, i + 1) = i*interp_coarse(end)^(i - 1);
    end  
end
[Q, R] = mgs(P_coarse);
Q = double(Q);

[coarse_Q_der, coarse_R_der] = mgs(P_coarse_der);
coarse_Q_tilde = (R*inv(coarse_R_der)*coarse_Q_der.').';
Q_tilde = real(double(coarse_Q_tilde));


% Forming Q_fine
P_fine = sym(zeros(n_over*(d - 1) + 1, d));
for i = 0 : d - 1
   P_fine(:, i + 1) = interp_fine.^i; 
end
Q_fine = P_fine*inv(R);

% Approximate the Gram polynomials by trigonometric polynomials using an
% SVD
fprintf( 'Performing SVD... \n'); 
if( mod(svd_order,2) == 0 )
    k_cos = sym(0 : svd_order/2 - 1);
    k_sin = sym(1 : svd_order/2 - 1);
else
    k_cos = sym(0 : (svd_order - 1)/2);
    k_sin = sym(1 : (svd_order - 1)/2);
end
X = [interp_fine; zero_fine];
C = [cos(2*sym('pi')*X*k_cos), sin(2*sym('pi')*X*k_sin)];
[U, S, V] = svd(C, 'econ');
Coeffs = sym(zeros(size(C, 2), d));
delta = diag(S);
delta_inv = 1./delta;


for i = 1 : d
    b = [Q_fine(:, i); sym(zeros(n_over*(Z - 1) + 1, 1))];
    Coeffs(:, i) = V*(delta_inv.*(U'*b)); % inverting the SVD
    r = max( abs( double( C*Coeffs(:, i) - b ) ) );
    fprintf( '\t%d\tresidual = %e\n', d, r );    
end

% Evaluating the trigonometric polynomials at the continuation points
fprintf( 'Evaluating at continuation points... \n '); 
A = double([cos(2*pi*continuation*k_cos), sin(2*pi*continuation*k_sin)]*Coeffs);

end