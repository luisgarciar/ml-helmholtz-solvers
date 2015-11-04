function f = rhs2(x,y)
% RHS2 Computes a discrete eigenfunction of the 
% Laplace operator on the square domain (0,1)^2
% to test a script for the Helmholtz equation
% INPUT: 
% [x,y]: meshgrid on (0,1)^2
% OUTPUT:
% f: eigenfunction evaluated on meshgrid x,y 
%%
m=5; n=5; k=0;
f = (sin(m*pi*x).*sin(n*pi*y));
f = ((n^2+m^2)*pi^2-k^2)*f;
end
