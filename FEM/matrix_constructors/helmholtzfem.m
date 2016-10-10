function [A] = helmholtzfem(k,np,eps)
%% HELMHOLTZFEM: Constructs a matrix for a 1D Helmholtz and Shifted Laplace model problem.
%  Constructs the matrix corresponding to the discretization of a 1D Helmholtz (or shifted Laplace)
%  problem
%       -u''- k^2*u-i*eps*u = f in (0,1)
%        u(0)       = 0   
%        u'(1)-ku(1)= 0
% 
% Use: A = helmholtzfem(k,np,eps)
%
% Note (For shifted Laplacian problems): 
% The imaginary shift is not multiplied by the term k^2
% as in the original publications of Erlangga et. al.
% For a convergent multigrid solver one needs eps~k^2 (eps=0.5k^2)
% For wave-independent GMRES convergence one needs eps~k 
% (but the MG solver will not converge)
%
% References:
% F. Ihlenburg and I. Babuska, Finite element solution of the 
% Helmholtz equation with high wave number Part I: The h-version of the FEM
%
% S Langdon , S Chandler-Wilde, Finite element methods for acoustic scattering 
% 
%  Input: 
%  k:   wavenumber
%  np:  size of the matrix (x=0 and interior points)
%  eps: imaginary shift    (eps=0 for Helmholtz problem)
%  
%  Output:
%  A:   FEM discretization matrix
%
%  Author: Luis Garcia Ramos, 
%          Institut fur Mathematik, TU Berlin
%  Version 0.1 - Jun 2016
%
%%%

%%


h  = 1/np;            %gridsize
e = ones(np,1);
A = (1/h)*spdiags([-e 2*e -e],[-1 0 1],np,np); A(np,np)=1/h; %Laplacian (Stiffness) Matrix
B = h*spdiags([(1/6)*e (2/3)*e (1/6)*e],[-1 0 1],np,np); B(np,np)=h/3; %Mass matrix
C = sparse(np,np); C(np,np)=1;  %Boundary term
A = A-(k^2+1i*eps)*B-1i*k*C;
end


