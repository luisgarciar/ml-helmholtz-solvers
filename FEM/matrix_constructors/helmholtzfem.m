function [A] = helmholtzfem(k,np,eps,bc)
%% HELMHOLTZFEM: Constructs a matrix for a 1D Helmholtz and Shifted Laplace model problem.
%  Constructs the matrix corresponding to the discretization of a 1D Helmholtz (or shifted Laplace)
%  problem with mixed or Sommerfeld boundary conditions
%
%       Mixed bc:
%       -u''- k^2*u-i*eps*u = f in (0,1)
%        u(0)       = 0   
%        u'(1)-ku(1)= 0
%
%       Sommerfeld bc:
%       -u''- k^2*u-i*eps*u = f in (0,1)
%        u'(0) + ku(0) = 0   
%        u'(1) - ku(1) = 0
%
%  Use: A = helmholtzfem(k,np,eps,bc)
%
%  Input: 
%  k:   wavenumber
%  np:  number of interior points in (0,1)
%  eps: imaginary shift    (eps=0 for Helmholtz problem)
%  bc: 'mix' or 'som'
%  
%  Output:
%  A:   FEM discretization matrix
%       size(A) = (np+2,np+2) if bc = 'som'  
%       size(A) = (np+1,np+1) if bc = 'mix'   
%
% Note (For shifted Laplacian problems): 
% The imaginary shift is not multiplied by the term k^2
% as in the original publications of Erlangga et. al.
% For a convergent multigrid solver one needs eps~k^2 (e.g. eps = 0.5k^2)
% For wave-independent GMRES convergence one needs eps~k 
% (but the MG solver will not converge)
%
% References:
% F. Ihlenburg and I. Babuska, Finite element solution of the 
% Helmholtz equation with high wave number Part I: The h-version of the FEM
%
% S Langdon, S Chandler-Wilde, Finite element methods for acoustic scattering 
% 
%  Author: Luis Garcia Ramos, 
%          Institut fur Mathematik, TU Berlin
%  Version 0.2 - Jun 2017
%%%%

switch bc
    case 'mix' %mixed boundary conditions
        h  = 1/(np+1);                 %gridsize
        e  = ones(np+1,1);
        A  = (1/h)*spdiags([-e 2*e -e],[-1 0 1],np+1,np+1); A(np+1,np+1)=1/h; %Laplacian Matrix
        B  = h*spdiags([(1/6)*e (2/3)*e (1/6)*e],[-1 0 1],np+1,np+1); B(np+1,np+1)=h/3; %Mass matrix
        C  = sparse(np+1,np+1); C(np+1,np+1)=1;  %Boundary term
        A  = A-(k^2+1i*eps)*B-1i*k*C;            %Helmholtz matrix

    case 'som' %sommerfeld boundary conditions
         h  = 1/(np+1);                 %gridsize
        e  = ones(np+2,1);
        A  = (1/h)*spdiags([-e 2*e -e],[-1 0 1],np+2,np+2); 
        A(1,1)=1/h; A(np+2,np+2)=1/h;       %Laplacian Matrix
        B  = h*spdiags([(1/6)*e (2/3)*e (1/6)*e],[-1 0 1],np+2,np+2); 
        B(1,1) = h/3; B(np+2,np+2)=h/3;     %Mass matrix
        C  = sparse(np+2,np+2); C(1,1)=1; C(np+2,np+2)=1;  %Boundary term
        A  = A-(k^2+1i*eps)*B-1i*k*C;            %Helmholtz matrix
end


