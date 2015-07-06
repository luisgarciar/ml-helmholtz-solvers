function f = eigflapl(x,y,m,n)
% EIGFLAPL Computes the (n,m) discrete eigenfunction of the 
%  Laplace operator on the square domain (0,1)^2
% INPUT: 
% [x,y]: meshgrid on (0,1)^2
% (n,m): index of the eigenfunction
% OUTPUT:
% f: eigenfunction evaluated on meshgrid x,y 
%%

f = sin(m*pi*x).*sin(n*pi*y);

end
