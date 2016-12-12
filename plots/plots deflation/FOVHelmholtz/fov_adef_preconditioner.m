%Field of values of finite element matrices
%ADEF preconditioner
clear all
close all

dim = 1;
k   = 500;
ppw = 15;
npc = 4;
eps = 0.5*k;
[npf,lev] = fem_npc_to_npf(npc,k,ppw);

A = helmholtzfem(k,npf,0); %Helmholtz matrix
M = helmholtzfem(k,npf,eps); %Shifted Laplace matrix

%eigv = eig(full(M\A));
%plot(real(eigv),imag(eigv),'.b','MarkerSize',22);
%axis([0 1 -0.5 0.5])
%xlabel('real(\lambda)','FontSize',14)
%ylabel('imag(\lambda)','FontSize',14)

R = fwrestrictionfem(npf,dim); 
Z = R'; %Prolongation. Deflation subspace: columns of Z
dim_def = size(Z,2);



%% Field of values of restricted matrix 
%To compute the field of values of the restricted matrix
%(ADEF) we need a basis Y for the orthogonal complement of
%the columns of Z, i.e., the nullspace of R.

Y = null(full(R));

%Restricted matrix
Ainvc = (Y'*(A\Y));
%Restricted ADEF matrix (only upper block)
B     = (Y'*(M\Y))*((Y'*(A\Y))\speye(length(Ainvc))); 
[FV_B, eigB] = fv(B,1,32);

%Explicit Deflation operator
Q = Z*((Z'*A*Z)\Z');
 
%Deflated-shifted Operator PadefA
PadefA =  M\(A-A*Z*((Z'*A*Z)\Z'*A))+ Z*((Z'*A*Z)\Z'*A);
 
%Field of values of PadefA
[FV_PadefA, eigPadefA] = fv(full(PadefA),1,32,1);

plot(real(FV_PadefA), imag(FV_PadefA),'b')     %Plot the field of Padef
axis('equal');
hold on
plot(real(eigPadefA), imag(eigPadefA), 'x')    %Plot the eigenvalues too.
plot(real(FV_B),imag(FV_B),'g')               %Plot the field of values of B
plot(real(eigB),imag(eigB), 'k.')    % Plot the eigenvalues too.
plot(0,0,'+k')

