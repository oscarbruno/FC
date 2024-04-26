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

function [Q, Q_tilde, A, ArQr, AlQl, ArQ_tilder, AlQ_tildel] =  ...
    precomp_fc_data(d, C, Z, E, n_over, modes_to_reduce, num_digits)


tol = 1e-16;
digits(num_digits);
N_coarse = sym(d + C + Z + E);
% modes_to_reduce = 4;
svd_order = sym(N_coarse - 2 * modes_to_reduce);



F1 = sym(zeros(d, d));
F2 = sym(zeros(C, C));

F1 = flipud(sym(eye(d)));
F2 = flipud(sym(eye(C)));

% Forming the grids we will use
N0 = sym(101);
h0 = sym(1/(N0 - 1));
x0 = sym(1 - (d - 1)*h0);
midpoint = sym(x0 + (d + C)*h0);


interp_coarse = sym(x0 + h0 * (0 : d - 1).'); 
interp_fine = sym(linspace(x0, 1, n_over*(d - 1) + 1).');
zero_fine = sym(midpoint + linspace(0, (Z - 1) * h0, n_over*(Z - 1) + 1).');
continuation = sym(x0 + h0 * (d : d - 1 + C).');


% Forming Q by using the modified Gram-Schmidt algorithm
P_coarse = sym(zeros(d, d));
P_coarse_der = sym(zeros(d, d));
for i = 0 : d - 1
    P_coarse(:, i + 1) = sym(interp_coarse.^i);
    P_coarse_der(1 : end - 1, i + 1) = sym(interp_coarse(1 : end - 1).^i);
    if (i == 0)
        P_coarse_der(end, i + 1) = 0;
    else
        P_coarse_der(end, i + 1) = sym(i*interp_coarse(end)^(i - 1));
    end  
end
[Q, R] = mgs(P_coarse);
% Q = double(Q);

[coarse_Q_der, coarse_R_der] = mgs(P_coarse_der);
coarse_Q_tilde = (R*inv(coarse_R_der)*coarse_Q_der.').';
Q_tilde = real((coarse_Q_tilde));


% Forming Q_fine
P_fine = sym(zeros(n_over*(d - 1) + 1, d));
for i = 0 : d - 1
   P_fine(:, i + 1) = sym(interp_fine.^i); 
end
Q_fine = (P_fine*inv(R));

% Approximate the Gram polynomials by trigonometric polynomials using an
% SVD
fprintf( 'Performing SVD... \n'); 
X = [interp_fine; zero_fine];
% X = [interp_fine; zero_fine] / ((N_coarse - 1) * h0);
% continuation = continuation / ((N_coarse - 1) * h0);


if( mod(svd_order,2) == 0 )
    k_cos = sym(0 : svd_order/2 - 1);
    k_sin = sym(1 : svd_order/2 - 1);
else
    k_cos = sym(0 : (svd_order - 1)/2);
    k_sin = sym(1 : (svd_order - 1)/2);
end

C = sym([cos(2*sym(pi)*X*k_cos / ((N_coarse - 1) * h0)), ...
        sin(2*sym(pi)*X*k_sin / ((N_coarse - 1) * h0))]);



[U, S, V] = svd(C, 'econ');
Coeffs = sym(zeros(size(C, 2), d));
delta = diag(S);
delta_inv = 1./delta;


for i = 1 : d
    b = [Q_fine(:, i); sym(zeros(n_over*(Z - 1) + 1, 1))];
    Coeffs(:, i) = V*(delta_inv.*(U'*b)); % inverting the SVD
    r = max( abs( vpa( C*Coeffs(:, i) - b ) ) );
    fprintf( '\t%d\tresidual = %e\n', i, r );    
end

% Evaluating the trigonometric polynomials at the continuation points
fprintf( 'Evaluating at continuation points... \n '); 


A = vpa([cos(2*sym(pi)*continuation * k_cos  / ((N_coarse - 1) * h0)),...
        sin(2*sym(pi)*continuation * k_sin  / ((N_coarse - 1) * h0))] * Coeffs);


ArQr = A*(Q.');
AlQl = F2*A*(Q.')*F1;
ArQ_tilder = A*(Q_tilde.');
AlQ_tildel = F2*A*(Q_tilde.')*F1;



end