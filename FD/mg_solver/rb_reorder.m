function P = rb_reorder(m)
%% RB_REORDERING Generates the permutation matrix P for red-black reordering 
% 
%   Use:    P = rb_reordering(m)  
%
%   Input: 
%       m:    number of interior points in 1-D grid
% 
%   Output: 
%        P:   mxm permutation matrix for red-black reordering     
%
% Input: 
%   m = number of interior grid points in each direction (m must be odd)
%
% Output: 
%   P = sparse permutation matrix 
% 
%   The vector x=(R1,B1,R2,B2,...) is reordered as
%   Px = (R1,R2,...,B1,B2,...).
%
%  The reordered matrix is A=P*A*P', and Ax=b
%  is transformed to (PAP')(Px)=Pb
%
%  Author: Luis Garcia Ramos, 
%          Institut fur Mathematik, TU Berlin
%          Version 0.1, Jun 2016

%%
if 2\(m+1) ~= round(2\(m+1)), 
  error('Input m must be odd.');
end

n = m^2;
P = speye(n,n); 

%setting up the permutation for even
%and odd indices
odd  = 1:1:(n+1)/2;
f    = @(x)(2*x-1); odd = f(odd);
even = 1:1:(n-1)/2;
g    = @(x) 2*x; even = g(even);

p = [odd even]; %permutation vector
P = P(p,:);     %permutation matrix

