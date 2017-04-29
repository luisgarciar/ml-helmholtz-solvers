function [eigvADEF1] = eigADEF1(k,np_int,eps)
%% EIGADEF1 Computes the eigenvalues of the 1D Helmholtz
%  operator with Dirichlet boundary conditions preconditioned by ADEF-1 
%  in combination with the shifted Laplacian
%
%   Input:
%       k: wavenumber
%       np_int: number of interior points
%       eps: shift
%
%   Output: 
%       eADEF1: nonzero eigenvalues of the Helmholtz matrix
%               preconditioned with ADEF-1
%
%
%   Author: Luis Garcia Ramos,          
%           Institut fur Mathematik, TU Berlin
%           Version 1.0, Mar 2017
%
% References: 
% Garcia R., Luis - Eigenvalue analysis of a two-Level deflation
% preconditioner for 1-D Helmholtz Problems 
% (unpublished notes, 2017)
%
%%
N = np_int;
if (mod(N+1,2)==1) 
    %print('invalid number of interior points')
    N = N+1; 
end

h = 1/(N+1);   %fine grid size 
n = (N+1)/2-1; %coarse grid size

j  = (1:N)';
g = pi*h*j;

%Notation
%A : Helmholtz matrix 
%M : shifted Laplacian, M = A-i*eps*I
%Pd: deflation operator 

%Eigenvalues of A
%eigvA = (2-2*cos(g1)-k^2*h^2)/h^2;

%Low frequency part of the spectrum of A
j      = (1:n)';
g      = (pi*h*j);
egvAj  = (2-2*cos(g)-k^2*h^2)/h^2;

%High frequency part of the spectrum of A
Nj  = N+1-j;
g   = (pi*h*Nj);
egvANj = (2-2*cos(g)-k^2*h^2)/h^2;

%Eigenvalues of M
%Low frequency part of the spectrum of M
g      = (pi*h*j);
egvMj  = (2-2*cos(g)-(k^2+1i*eps)*h^2)/h^2;

%High frequency part of the spectrum of M
Nj  = N+1-j;
g   = (pi*h*Nj);
egvMNj = (2-2*cos(g)-(k^2+1i*eps)*h^2)/h^2;

% We compute the eigenvalues of M.Pd.A
% Each eigenvalue is of the form xj = yj*zj  
% yj: weighted Harmonic mean of eigenvalues of inv(A)
% zj: convex combination of eigenvalues of inv(M)

cj = cos(pi*h*j/2);  sj = sin(pi*h*j/2);
aj = (cj.^4)./(cj.^4+sj.^4);  bj = (sj.^4)./(cj.^4+sj.^4);

%Eigenvalues of Helmholtz matrix preconditioned by ADEF-1
yj = 1./(aj./egvANj + bj./egvAj);
zj = (aj./egvMNj + bj./egvMj);

eigvADEF1 = yj.*zj;

end

