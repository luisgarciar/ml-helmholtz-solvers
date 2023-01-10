%Eigenvalues of 1-D Helmholtz matrices preconditioned by
%ADEF-1 and the shifted Laplacian%
clear all; close all; clc
k   = 50;
ppw = 40;
eps = 0.5*k^2;
b1  = 1; b2=0.5;
npf = ceil(ppw*k/(2*pi))-1;
dim = 1;
bc = 'dir';

if (mod(npf+1,2)==1) 
    npf = npf+1; 
end
npc = (npf-1)/2;

% Eigenvalues of ADEF-1 computed symbolically using the new formulas
[eig_adef1] = eigADEF1(k,npf,b1,b2); 

%Helmholtz and shifted Laplacian Matrices
A = helmholtz(k,0,npf,bc);
M = helmholtz(k,eps,npf,bc);
S = M\A;
eigS = eig(full(S)); 
eigS = sort(eigS); %Eigenvalues of CSL

%Intergrid operators and deflation (ADEF-1) preconditioner
Y = fwrestriction(npf,dim,bc);
Z = lininterpol(npc,dim,bc);
E = Z'*A*Z;   %coarse grid operator
[L,U] = lu(E); I = eye(length(A));
P = M\(I-A*Z*(E\Z'))+Z*(E\Z');
PadefA  = P*A;  %deflated operator;

%Eigenvalues of linear system preconditioned by ADEF-1 
%computed numerically (full linear system)
[FV_ADEF1, eigvADEF1] = fv(full(PadefA),1,32,1);

%Eigenvalues computed numerically using the restricted system
V = null(full(Z')); %basis vectors for the orthogonal complement of range of Z
B = (V'*(A\V))\(V'*(M\V)); %restricted matrix
[FV_B, eigvB] = fv(full(B),1,40,1);


%Plots
figure(1);
plot(real(eig_adef1),imag(eig_adef1),'*b','MarkerSize',12);
axis([-1 2 -1 2]); %grid on
title('Eigenvalues computed from the formulas')

figure(2);
plot(real(eigvB),imag(eigvB),'*r','MarkerSize',12);
title('Eigenvalues and FOV computed from the restricted matrix')
axis([-1 2 -1 2]); %grid on
hold on
plot(real(FV_B),imag(FV_B),'b') %Plot the field of values of B


figure(3)
plot(real(eigvADEF1),imag(eigvADEF1),'+k','MarkerSize',12);
title('Eigenvalues and FOV computed from the full matrix')
%axis([-1 2 -1 2]); %grid on
hold on
plot(real(FV_ADEF1),imag(FV_ADEF1),'b') %Plot the field of values of B


