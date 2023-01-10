%test of complex gmress with custom givens rotations
%
n = 10;
A = randn(n,n);
b = randn(n,1);
M = eye(n);
x0 = zeros(n,1);
restrt = 20;
maxit = 20;
tol   = 1e-10;

[x, error, iter, flag] = gmress(A,x0,b,M,restrt,maxit,tol);

relres = norm(b-A*x)/norm(b)

error

