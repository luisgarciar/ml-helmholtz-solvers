%% csminres_helmholtz
%% In this script we compare the performance of the (unpreconditioned) GMRES
%% and CS-MINRES methods for various problems.

%% Random b
n = 30;
b = rand(n,1);
b = sort(b,'descend');

%% Constant f
b = ones(n,1);

%% Eigenvalues of the matrix
lambda = rand(n,1);
lambda = cos(lambda)+1i*sin(lambda);

%% Construction of unitary, complex symmetric matrix A
[U,~] = qr(randn(n,n));
    A = U*diag(lambda)*U';
    

%% Running GMRES and CSMINRES    
[X,FLAG,RELRES,ITER,RESVEC]  = gmres(A,b,[],1e-12,n);
[x,flag,iter,~,~,relres,~,~,~,~,~,resvec,aresvec] = csminresqlp(A,b);

semilogy([0:length(RESVEC)-1],RESVEC/norm(b),'b-.',[0:length(resvec)-1],resvec/norm(b),'r-+');
%semilogy([0:length(resvec)-1],resvec/norm(b),'r-+');

legend('Residuals')
xlabel('Iteration')
ylabel('Relative residuals')

 