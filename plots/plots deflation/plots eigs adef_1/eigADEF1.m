function [eigvADEF1] = eigADEF1(k,ppw,b1,b2)
%% EIGADEF1 Generates the eigenvalues of the 1D Helmholtz
%  operator with Dirichlet boundary conditions
%  preconditioned by the ADEF-1 preconditioner
%  in combination with the shifted Laplacian
%
%   Input:
%       k: Wavenumber
%       ppw: points per wavelength
%       (b1,b2): shift
%
%   Output: 
%           eADEF1: nonzero eigenvalues of the preconditioned system           
%
%   Author: Luis Garcia Ramos,          
%           Institut fur Mathematik, TU Berlin
%           Version 0.1, May 2016
%
%%
N = ceil(ppw*k/(2*pi))-1;
if (mod(N+1,2)==1) 
    N = N+1; 
end

h = 1/(N+1);   %meshgrid size 
n = (N+1)/2-1; %Coarse grid size

j  = (1:N)';
g1 = pi*h*j;

%Notation
%A : Helmholtz matrix 
%M : shifted Laplacian
%Pd: deflation operator 

%Eigenvalues of A
%eigvA = (2-2*cos(g1)-k^2*h^2)/h^2;

%Low frequency part of the spectrum of A
j   = (1:n)';
g   = (pi*h*j);
egvAj  = (2-2*cos(g)-k^2*h^2)/h^2;

%High frequency part of the spectrum of A
Nj  = N+1-j;
g   = (pi*h*Nj);
egvANj = (2-2*cos(g)-k^2*h^2)/h^2;

%Eigenvalues of M
%eigvM = (2-2*cos(g1)-k^2*(b1-1i*b2)*h^2)/h^2;

%Low frequency part of the spectrum of M
g      = (pi*h*j);
egvMj  = (2-2*cos(g)-k^2*(b1-1i*b2)*h^2)/h^2;
%lj = (sin(g).^2-(k*h)^2)./(sin(g).^2-(k*h)^2*(b1-1i*b2));

%High frequency part of the spectrum of M
Nj  = N+1-j;
g   = (pi*h*Nj);
egvMNj = (2-2*cos(g)-k^2*(b1-1i*b2)*h^2)/h^2;

% We compute the eigenvalues of M.Pd.A
% Each eigenvalue is of the form xj = yj.zj  
% yj: weighted Harmonic mean of eigenvalues of inv(A)
% zj: convex combination of eigenvalues of inv(M)

cj = cos(pi*h*j/2);  sj = sin(pi*h*j/2);
aj = (cj.^4)./(cj.^4+sj.^4);  bj = (sj.^4)./(cj.^4+sj.^4);

%Eigenvalues of Helmholtz preconditioned by ADEF-1
yj = 1./(aj./egvANj + bj./egvAj);
zj = (aj./egvMNj + bj./egvMj);

eigvADEF1 = yj.*zj;

end

