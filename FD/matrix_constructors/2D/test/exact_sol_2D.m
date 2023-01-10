function [u] = exact_sol_2D(x,y,m,n)
%Computes an exact solution to the 2D Helmholtz problem
% u''-k^2u = sin(m*pi*x)*sin(n*pi*y) in (0,1)x(0,1)
% u = 0 on the boundary
%
% Use: [u] = exact_sol_2D(x,y,k,m,n)
%
% Input:  x,y:  arrays in (0,1)x(0,1)
%         k,m,n: Wavenumber and parameters
%
% Output: u = f(x,y);
% 
%
assert(size(x)==size(y),'x,y must be the same size')
u = sin(m*pi*x).*sin(n*pi*y);


end

