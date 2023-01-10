function [eCSL,eDCSL] = eigSL(k,ppw,b1,b2)
%% EIGSL Generates the eigenvalues of the 1D Helmholtz
%  operator with Dirichlet boundary conditions
%  preconditioned by the shifted Laplacian and
%  the deflated CSL 
%
%   Input:
%       k: Wavenumber
%       ppw: points per wavelength
%       (b1,b2): shift
%
%   Output: 
%           eCSL
%
%           eDCSL
%
%
%   Author: Luis Garcia Ramos,          
%           Institut fur Mathematik, TU Berlin
%           Version 0.1, Feb 2016
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

%Eigenvalues of shifted Laplacian
eCSL  = (2-2*cos(g1)-k^2*h^2)./(2-2*cos(g1)-k^2*h^2*(b1-1i*b2));

%Smooth part of the spectrum of shifted Laplacian
j   = (1:n)';
g   = (pi*h*j);
lj  = (2-2*cos(g)-k^2*h^2)./(2-2*cos(g)-k^2*h^2*(b1-1i*b2));
%lj = (sin(g).^2-(k*h)^2)./(sin(g).^2-(k*h)^2*(b1-1i*b2));

%Oscillatory part of the spectrum of shifted Laplacian
Nj  = N+1-j;
g   = (pi*h*Nj);
lNj = (2-2*cos(g)-k^2*h^2)./(2-2*cos(g)-k^2*h^2*(b1-1i*b2));

%Computing the eigenvalues of the deflated shifted Laplacian
% Each eigenvalue of the deflated CSL is a weighted Harmonic mean 
%of a smooth and an oscillatory eigenvalue of the CSL
cj = cos(pi*h*j/2);      sj = sin(pi*h*j/2);
aj = (cj.^4)./(cj.^4+sj.^4);  bj = (sj.^4)./(cj.^4+sj.^4);

%Eigenvalues of deflated CSL
eDCSL = 1./(aj./lNj+bj./lj);

end

