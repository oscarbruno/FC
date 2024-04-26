
% A Modified Gram-Schmidt routine to compute the QR decomposition.
%
% Input
%   A : m x n sized matrix
%   
% Outputs:
%   Q : m x n matrix of the QR decomposition of A
%   R : n x n matrix of the QR decomposition of A
%

function [Q, R] = mgs(A)
m = size(A,1);
n = size(A,2);
Q = sym(zeros(m,n));
R = sym(zeros(n,n));

% R(1,1) = norm(A(:,1));
R(1,1) = sqrt(A(:,1).' * A(:,1));

Q(:,1) = A(:,1)/R(1,1);

% for j=2:n
%     q = A(:,j);
%     for i=1:j-1
%         R(i,j) = q'*Q(:,i);
%         q = q - R(i,j)*Q(:,i);
%     end
% 
% %     R(j,j) = norm(q);
%     R(j, j) = sqrt(q.' * q);
%     if (R(j,j)==0)
%         break;
%     else
%         Q(:,j) = q/R(j,j);
%     end
% end
% end

    for j=2:n
        Q(:, j) = A(:,j);
        for i=1:j-1
            R(i,j) = Q(:, j).'*Q(:,i);
            Q(:, j) = Q(:, j) - R(i,j)*Q(:,i);
        end
        % (Full) re-orthogonalization; do the above procedure again to
        % safuguard against roundoff errors and ensure orthoganality
        for i = 1 : j-1
            proj_subspace = Q(:,i).' * Q(:,j); 
            R(i, j) = R(i, j) + proj_subspace;        
            Q(:,j) = Q(:,j) - proj_subspace * Q(:,i);
        end 
    %     R(j,j) = norm(q);
        R(j, j) = sqrt(Q(:, j).' * Q(:, j));
        Q(:,j) = Q(:, j)/R(j,j);
    %     if (R(j,j)==0)
    %         break;
    %     else
    %         Q(:,j) = Q(:, j)/R(j,j);
    %     end
    end
end

% for ii = 1:d
%     Q(:,ii) = std_basis(:,ii);    
%     for jj = 1:ii-1
%         proj_subspace = Q(:,jj).' * Q(:,ii);        
%         R(jj, ii) = proj_subspace;       
%         Q(:,ii) = Q(:,ii) - proj_subspace * Q(:,jj);
%     end    
%     % (Full) re-orthogonalization; do the above procedure again to
%     % safuguard against roundoff errors and ensure orthoganality
%     for jj = 1:ii-1
%         proj_subspace = Q(:,jj).' * Q(:,ii); 
%         R(jj, ii) = R(jj, ii) + proj_subspace;        
%         Q(:,ii) = Q(:,ii) - proj_subspace * Q(:,jj);
%     end    
%     % normalize for orthonormality
%     p_nrm = sqrt(Q(:,ii).' * Q(:,ii));
%     Q(:,ii) = Q(:,ii) / p_nrm;
%     R(ii, ii) = p_nrm;
% end
    
  