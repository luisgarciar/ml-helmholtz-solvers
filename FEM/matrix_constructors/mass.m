function [M] = mass(np)
%% MASS: Constructs a mass matrix for the 1D Helmholtz model problem.
%  Constructs the mass corresponding to the discretization of a problem
%  with P1 finite elements in (0,1). The model problem is given by
%
%       -u''- k^2*u-i*eps*u = f in (0,1)
%        u(0)       = 0   
%        u'(1)-ku(1)= 0
% 
% Use: A = mass(np)
%
% References:
% F. Ihlenburg and I. Babuska, Finite element solution of the 
% Helmholtz equation with high wave number Part I: The h-version of the FEM
%
% S Langdon , S Chandler-Wilde, Finite element methods for acoustic scattering 
% 
%  Input: 
%  np:  size of the matrix (x=0 and interior points)
%  
%  Output:
%  M:   Mass matrix
%
%  Author: Luis Garcia Ramos, 
%          Institut fur Mathematik, TU Berlin
%  Version 0.1 - Jun 2016
%

%%
h = 1/np;            %gridsize
e = ones(np,1);
M = h*spdiags([(1/6)*e (2/3)*e (1/6)*e],[-1 0 1],np,np); M(np,np)=h/3; %Mass matrix
end


