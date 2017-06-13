%Field of values of finite element matrices 1D
%ADEF preconditioner
clear all
close all

dim = 1;
k   = 20;
ppw = 0.1;
npc = 4;
eps = 100000*k^2;
[npf,lev] = fem_npc_to_npf(npc,k,ppw);

A  = helmholtzfem(k,npf,0); %Helmholtz matrix
SL = helmholtzfem(k,npf,eps); %Shifted Laplace matrix

Ahat = full(SL\A);
eigv = eig(Ahat);
plot(real(eigv),imag(eigv),'.b','MarkerSize',22);
axis([0 1 -0.5 0.5])
xlabel('real(\lambda)','FontSize',14)
ylabel('imag(\lambda)','FontSize',14)

R = fwrestrictionfem(npf,dim);
Z = R'; %Prolongation. Deflation subspace: columns of Z
dim_def = size(Z,2);

%Deflated-shifted Operator PadefA
M = mass(npf);
sqrtM = sqrtm(full(M));
I = eye(length(A));
P = full(SL\(I-A*Z*((Z'*A*Z)\Z')));
P = sqrtM*P*sqrtM;

%Field of values of Ahat Sl-preconditioned Helmholtz matrix
[fvAhat, eigAhat] = fv(Ahat,1,32,1);

%Field of values of PadefA in Minv inner product
[fvAP] = 1+1i*eps*fv(P,1,32,1);

%Plots
figure(2)
plot(real(fvAhat), imag(fvAhat),'b') %Plot the FOV of Padef
hold on
plot(real(eigAhat), imag(eigAhat),'r+')    %Plot the eigenvalues too.
axis('equal');
plot(0,0,'+k')
plot(real(fvAP), imag(fvAP),'k*') %Plot the FOV of Padef


%hold on
%plot(real(eigAP), imag(eigAP),'r+')    %Plot the eigenvalues too.

% 
% %% 1D Dirichlet Example
% dim = 1;
% npc = 4;
% bc = 'dir';
% %wavenumber and imaginary shift of shifted Laplacian
% k   = 200;  eps = 0.5*k^2; %Helmholtz problem
% ppw = 15;   %number of points per wavelength
% [npf,numlev] = fd_npc_to_npf(npc,k,ppw);  %number of points in finest grid (1D)
% 
% A = helmholtz(k,0,npf,bc); 
% M = helmholtz(k,eps,npf,bc);  %Shifted Laplace matrix
% 
% R = fwrestriction(npf,dim,bc); 
% Z = R'; %Prolongation. Deflation subspace: columns of Z
% dim_def = size(Z,2);
% 
% % Field of values of restricted matrix 
% %To compute the field of values of the restricted matrix
% %(ADEF) we need a basis Y for the orthogonal complement of
% %the columns of Z, i.e., the nullspace of R.
% 
% Y = null(full(R));
% 
% %Restricted matrix
% Ainvc = (Y'*(A\Y));
% %Restricted ADEF matrix (only upper block)
% B     = (Y'*(M\Y))*((Y'*(A\Y))\speye(length(Ainvc))); 
% [FV_B, eigB] = fv(B,1,32);
% 
% %Explicit Deflation operator
% Q = Z*((Z'*A*Z)\Z');
% %Deflated-shifted Operator PadefA
% PadefA =  M\(A-A*Z*((Z'*A*Z)\Z'*A))+ Z*((Z'*A*Z)\Z'*A);
% 
% [U,R] = qr(Z,0);
% Porth = U*U';
% 
% I = eye(length(A));
% %PadefAcorr = PadefA-Porth*PadefA*(I-Porth);
% PadefAcorr = (I-Porth)*PadefA + Porth;
% 
% %Plots
% %Field of values of PadefA
% [FV_PadefA, eigPadefA] = fv(full(PadefA),1,32,1);
% [FV_PadefAcorr, eigPadefAcorr] = fv(full(PadefAcorr),1,32,1);
% 
% 
% plot(real(FV_PadefA), imag(FV_PadefA),'b') %Plot the FOV of Padef
% %plot(real(eigPadefA), imag(eigPadefA), 'r+')    %Plot the eigenvalues too.
% axis('equal');
% hold on
% %plot(real(FV_PadefAcorr), imag(FV_PadefAcorr),'k')   %Plot the FOV of Padefcorr
% %plot(real(eigPadefAcorr), imag(eigPadefAcorr), 'kx') %Plot the eigenvalues too.
% plot(real(FV_B),imag(FV_B),'r')     %Plot the field of values of B
% %plot(real(eigB),imag(eigB), 'k.')  %Plot the eigenvalues too.
% %plot(0,0,'+k')

% 
% %% 2D Sommerfeld Example
% clear all
% dim = 2;
% npc = 4;
% bc = 'som';
% %wavenumber and imaginary shift of shifted Laplacian
% k   = 10;  eps = 0.5*k^2; %Helmholtz problem
% ppw = 10;   %number of points per wavelength
% [npf,numlev] = fd_npc_to_npf(npc,k,ppw);  %number of points in finest grid (1D)
% 
% A = helmholtz2(k,0,npf,npf,bc); 
% M = helmholtz2(k,eps,npf,npf,bc);  %Shifted Laplace matrix
% 
% R = fwrestriction(npf,dim,bc); 
% Z = R'; %Prolongation. Deflation subspace: columns of Z
% dim_def = size(Z,2);
% 
% % Field of values of restricted matrix 
% %To compute the field of values of the restricted matrix
% %(ADEF) we need a basis Y for the orthogonal complement of
% %the columns of Z, i.e., the nullspace of R.
% 
% Y = null(full(R));
% 
% %Restricted matrix
% Ainvc = (Y'*(A\Y));
% %Restricted ADEF matrix (only upper block)
% B     = (Y'*(M\Y))*((Y'*(A\Y))\speye(length(Ainvc))); 
% [FV_B, eigB] = fv(B,1,32);
% 
% %Explicit Deflation operator
% Q = Z*((Z'*A*Z)\Z');
% %Deflated-shifted Operator PadefA
% PadefA =  M\(A-A*Z*((Z'*A*Z)\Z'*A))+ Z*((Z'*A*Z)\Z'*A);
% 
% [U,R] = qr(Z,0);
% Porth = U*U';
% 
% I = eye(length(A));
% %PadefAcorr = PadefA-Porth*PadefA*(I-Porth);
% PadefAcorr = (I-Porth)*PadefA + Porth;
% 
% %Plots
% %Field of values of PadefA
% [FV_PadefA, eigPadefA] = fv(full(PadefA),1,32,1);
% [FV_PadefAcorr, eigPadefAcorr] = fv(full(PadefAcorr),1,32,1);
% 
% 
% plot(real(FV_PadefA), imag(FV_PadefA),'b') %Plot the FOV of Padef
% hold on
% plot(real(eigPadefA), imag(eigPadefA), 'r+')    %Plot the eigenvalues too.
% axis('equal');
% plot(real(FV_PadefAcorr), imag(FV_PadefAcorr),'k')   %Plot the FOV of Padefcorr
% plot(real(eigPadefAcorr), imag(eigPadefAcorr), 'kx') %Plot the eigenvalues too.
% plot(real(FV_B),imag(FV_B),'gx')     %Plot the field of values of B
% plot(real(eigB),imag(eigB), 'k.')  %Plot the eigenvalues too.
% plot(0,0,'+k')
% 
% 
% %% Testing if there is a difference between the matrices in GMRES
% tol   = 1e-10;
% maxit = length(A);
% b     = randn(length(A),1);
% 
% [x1,flag1,relres1,iter1,resvec1] = gmres(PadefA,b,[],tol,maxit,[]);
% [x2,flag2,relres2,iter2,resvec2] = gmres(PadefAcorr,b,[],tol,maxit,[]);
% 
% iter1
% iter2
% 
% figure(3)
% semilogy(1:(iter2(2)+1),resvec2'/resvec2(1),'r-+')
% hold on
% semilogy(1:(iter1(2)+1),resvec1'/resvec1(1),'b-+')