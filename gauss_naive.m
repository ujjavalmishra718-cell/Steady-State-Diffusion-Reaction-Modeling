function x = gauss_naive(A,b)
% x = GaussNaive(A,b): Gauss elimination without pivoting.
% input:
% A = coefficient matrix
% b = right hand side vector (non-homogeneity)
% output:
% x = solution vector
[m,n] = size(A);
if m~=n, error('Matrix A must be square'); end
nb = n+1;
B = [A b];      % Augmented matrix
% forward elimination
for k = 1:n-1
    for i = k+1:n
        factor = B(i,k)/B(k,k);
        B(i,k:nb) = B(i,k:nb)-factor*B(k,k:nb);
    end
end
% back substitution
x = zeros(n,1);
x(n) = B(n,nb)/B(n,n);
for i = n-1:-1:1
    x(i) = (B(i,nb)-B(i,i+1:n)*x(i+1:n))/B(i,i);
end
end