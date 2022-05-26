
% A routine used to pre-compute the matrix products necessary for the
% building of the Fourier Continuation


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