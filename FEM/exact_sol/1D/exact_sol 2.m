function [u] = exact_sol(k,x)
%EXACT_SOL  Computes the solution of the 1D Helmholtz 
%           problem
%       -u''- k^2*u = 1 in (0,1)
%        u(0)       = 0   
%        u'(1)-ku(1)= 0
% 
%Use: u = exact_sol(k,x)
%
% Input: 
%  k:   wavenumber 
%  x:   vector in (0,1)
%
% Output:
%  u:  value of u(x)
%
% References:
% F. Ihlenburg and I. Babuska, Finite element solution of the 
% Helmholtz equation with high wave number Part I: The h-version of the FEM 
%
% Author: Luis Garcia Ramos, 
%          Institut fur Mathematik, TU Berlin
% Version 0.1 - Jun 2016
%
%
u = (-1i*exp(1i*k).*sin(k*x)+exp(1i*k*x).*(1i*sin(k*x)-cos(k*x)+1))./k^2;
end

