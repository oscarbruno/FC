
% A routine used to build save the matrices necessary for the
% generation of the Fourier Continuation.
%
%
% Author:
%
%   Daniel Leibovici
%   email: dleibovi@caltech.edu
%

function generate_bdry_continuations(d, C, E, Z, n_ovr, num_digits)

tic;
fprintf('Performing precomputations...\n');
[Q, Q_tilde, A] =  precomp_fc_data(d, C, Z, E, n_ovr, num_digits);
save(['FC_data/A_d',num2str(d),'_C_', num2str(C), '.mat'], 'A');
save(['FC_data/Q_d',num2str(d),'_C_', num2str(C), '.mat'], 'Q');
save(['FC_data/Q_tilde_d',num2str(d),'_C_', num2str(C), '.mat'], 'Q_tilde');
toc;

end
