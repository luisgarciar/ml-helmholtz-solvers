function f = dirac1d(x)
% DIRAC1D Computes a right hand side vector with a 
% Dirac delta function in 1D
%
%  USAGE: f = rhs2d(x,h)
%
%  INPUT: 
%  x: grid on (0,1)
%   
%  OUTPUT:
%  f: discrete Dirac function centered on the midpoint of the grid
%
%%

f = zeros(size(x));
h = x(2)-x(1);
f(ceil(length(f)/2)) = 1/h;
        
    
end
