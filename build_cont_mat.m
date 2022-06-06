
% A routine used to pre-compute the matrix products necessary for the
% building of the Fourier Continuation
%
% Inputs: 
%   A : matrix containing the values of the extended Gram polynomials at 
%   the continuation points
%   Q : matrix containing the values of the Gram polynomials (for Dirichlet
%   boundary conditions)at the d points of the interval. 
%   Q_tilde : matrix containing the values of the Gram polynomials (for
%   Neumann boundary conditions)at the d points of the interval. 
%
% Outputs:
%   ArQr : matrix A_r*Q_r^T in eq (25) in [1]
%   AlQl : matrix A_l*Q_l^T in eq (25) in [1]
%   ArQ_tilder : matrix A_r*\tilde{Q}_r^T in eq (26) in [1]
%   AlQ_tildel : matrix A_r*\tilde{Q}_r^T in eq (26) in [1]
%
% Author:
%   Daniel Leibovici
%   email: dleibovi@caltech.edu
%


function [ArQr, AlQl, ArQ_tilder, AlQ_tildel] = build_cont_mat(A, Q, Q_tilde)

d = size(Q_tilde, 1);
C = size(A, 1);
F1 = flipud(eye(d));
F2 = flipud(eye(C));
ArQr = A*(Q.');
AlQl = F2*A*(Q.')*F1;
ArQ_tilder = A*(Q_tilde.');
AlQ_tildel = F2*A*(Q_tilde.')*F1;


end