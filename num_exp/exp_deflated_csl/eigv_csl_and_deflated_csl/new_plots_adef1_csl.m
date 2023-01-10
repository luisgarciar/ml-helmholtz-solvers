%Eigenvalues of 1-D Helmholtz matrices preconditioned by
%ADEF-1 and the shifted Laplacian%
clear all; close all; clc
k   = 100;
eps = 0.5*k^2;
b1  = 1; b2=1;
npf = ceil(ppw*k/(2*pi))-1;
dim = 1;
bc = 'dir';

if (mod(npf+1,2)==1) 
    npf = npf+1; 
end
npc = (npf-1)/2;

% Eigenvalues of ADEF-1 computed symbolically using the new formulas
[eig_adef1] = eigADEF1(k,ppw,b1,b2); 
figure(1);
plot(real(eig_adef1),imag(eig_adef1),'*b','MarkerSize',12);
axis([-1 2 -1 0.5]); %grid on
hold on

%% Eigenvalues of ADEF-1 computed numerically constructing the matrices
A = helmholtz(k,0,npf,bc);
M = helmholtz(k,eps,npf,bc);
S = M\A;
eigS = eig(full(S)); 
eigS = sort(eigS); %Eigenvalues of CSL
figure(2)
plot(real(eigS),imag(eigS),'*b','MarkerSize',12);


Y = fwrestriction(npf,dim,bc);
Z = lininterpol(npc,dim,bc);
E = Z'*A*Z;   %coarse grid operator
[L,U] = lu(E); I = eye(length(A));
P = I-A*Z*(E\Z')+Z*(E\Z');
PA = P*A;  %deflated operator;
MinvPA  = M\(PA); %preconditioned deflated operator

%Eigenvalues of linear system preconditioned by ADEF-1 computed numerically
eigvADEF1 = eig(MinvPA); 
eigvADEF1 = sort(eigvADEF1,'descend'); %
%eigvADEF1 = eigvADEF1(1:length(eigvADEF1),1); %Nonzero eigenvalues

figure(3)
plot(real(eigvADEF1),imag(eigvADEF1),'+r','MarkerSize',12);
