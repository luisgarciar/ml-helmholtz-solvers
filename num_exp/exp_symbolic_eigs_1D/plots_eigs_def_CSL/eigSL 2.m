function [eCSL,eDCSL] = eigSL(k,N,b1,b2)
%% EIGSL Generates the eigenvalues of the 1D Helmholtz operator 
% with Dirichlet boundary conditions  preconditioned by the shifted Laplacian 
% and in combination with deflation
%
%   Input:
%       k: Wavenumber
%       N: number of interior grid points (size of linear system)
%      (b1,b2): shift
%
%   Output: 
%           eCSL : eigenvalues of shifted Laplacian
%           eDCSL: eigenvalues of deflated shifted Laplacian
%
%   Author: Luis Garcia Ramos,          
%           Institut fur Mathematik, TU Berlin
%           Version 1.0, Dec 2016
%
%Reference:
%
%
%%

if (mod(N+1,2)==1) 
    N = N+1; 
end

h = 1/(N+1);   %meshgrid size 
n = (N+1)/2-1; %Coarse grid size

j  = (1:N)';
g1 = pi*h*j;

%Eigenvalues of shifted Laplacian
eCSL  = (2-2*cos(g1)-k^2*h^2)./(2-2*cos(g1)-k^2*h^2*(b1+1i*b2));

%Smooth part of the spectrum of shifted Laplacian
j   = (1:n)';
g   = (pi*h*j);
lj  = (2-2*cos(g)-k^2*h^2)./(2-2*cos(g)-k^2*h^2*(b1+1i*b2));
%lj = (sin(g).^2-(k*h)^2)./(sin(g).^2-(k*h)^2*(b1-1i*b2));

%Oscillatory part of the spectrum of shifted Laplacian
Nj  = N+1-j;
g   = (pi*h*Nj);
lNj = (2-2*cos(g)-k^2*h^2)./(2-2*cos(g)-k^2*h^2*(b1+1i*b2));

%Computing the eigenvalues of the deflated shifted Laplacian
%Each eigenvalue of the deflated CSL is a weighted Harmonic mean 
%of a smooth and an oscillatory eigenvalue of the CSL
cj = cos(pi*h*j/2);           sj = sin(pi*h*j/2);
aj = (cj.^4)./(cj.^4+sj.^4);  bj = (sj.^4)./(cj.^4+sj.^4);

%Eigenvalues of deflated CSL
eDCSL = 1./(aj./lNj+bj./lj);

end

