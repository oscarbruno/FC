
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

R(1,1) = norm(A(:,1));
Q(:,1) = A(:,1)/R(1,1);

for j=2:n
    q = A(:,j);
    for i=1:j-1
        R(i,j) = q'*Q(:,i);
        q = q - R(i,j)*Q(:,i);
    end
    R(j,j) = norm(q);
    if (R(j,j)==0)
        break;
    else
        Q(:,j) = q/R(j,j);
    end
end
end