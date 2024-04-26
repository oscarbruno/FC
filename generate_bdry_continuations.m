
% A routine used to build save the matrices necessary for the
% generation of the Fourier Continuation.
%
%
% Author:
%
%   Daniel Leibovici
%   email: dleibovi@caltech.edu
%

function generate_bdry_continuations(d, C, E, Z, n_ovr, modes_to_reduce, ...
    num_digits)

tic;
fprintf('Performing precomputations...\n');
[Q, Q_tilde, A, ArQr, AlQl, ArQ_tilder, AlQ_tildel] =  ... 
    precomp_fc_data(d, C, Z, E, n_ovr, modes_to_reduce, num_digits);

save(['FC_data/A_d',num2str(d),'_C', num2str(C), '.mat'], 'A');
save(['FC_data/Q_d',num2str(d),'_C', num2str(C), '.mat'], 'Q');
save(['FC_data/Q_tilde_d',num2str(d),'_C', num2str(C), '.mat'], 'Q_tilde');


% save(['FC_data/ArQr_d',num2str(d),'_C', num2str(C), '.mat'], 'ArQr');
% save(['FC_data/AlQl_d',num2str(d),'_C', num2str(C), '.mat'], 'AlQl');
% save(['FC_data/ArQ_tilder_d',num2str(d),'_C', num2str(C), '.mat'], 'ArQ_tilder');
% save(['FC_data/AlQ_tildel_d',num2str(d),'_C', num2str(C), '.mat'], 'AlQ_tildel');
% 
% 
% dlmwrite(['FC_data/A', num2str(d), 'C', num2str(C), '.dat'], A, ...
%          'delimiter', ',', 'precision', 30);
% dlmwrite(['FC_data/Q', num2str(d), 'C', num2str(C), '.dat'], Q, ...
%          'delimiter', ',', 'precision', 30);
% dlmwrite(['FC_data/Q_tilde', num2str(d), 'C', num2str(C), '.dat'], Q_tilde, ...
%  'delimiter', ',', 'precision', 30);
% dlmwrite(['FC_data/ArQr', num2str(d), 'C', num2str(C), '.dat'], ArQr, ...
%      'delimiter', ',', 'precision', 30);
% dlmwrite(['FC_data/AlQl', num2str(d), 'C', num2str(C), '.dat'], AlQl, ...
%      'delimiter', ',', 'precision', 30);   
% dlmwrite(['FC_data/ArQ_tilder', num2str(d), 'C', num2str(C), '.dat'], ArQ_tilder, ...
%      'delimiter', ',', 'precision', 30); 
% dlmwrite(['FC_data/AlQ_tildel', num2str(d), 'C', num2str(C), '.dat'], AlQ_tildel, ...
%      'delimiter', ',', 'precision', 30);       
    
toc;

end
