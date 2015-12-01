function f=rhs2d(x,y)
%% RHS2D Computes a discrete eigenfunction of the 
% Laplace operator on the square domain (0,1)^2
% to test a script for the Helmholtz equation
%
% Use: f = rhs2d(x,y)
%
% Input: 
% [x,y]: meshgrid on (0,1)^2
%
% Output:
% f: eigenfunction evaluated on meshgrid x,y 
%%
m=5; n=5; k=0;
f = (sin(m*pi*x).*sin(n*pi*y));
f = ((n^2+m^2)*pi^2-k^2)*f;
end
