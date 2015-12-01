function k = klayer(x,y,kref)
%% KLAYER Constructs a matrix of layered wavenumbers for the Helmholtz problem
%        
% Use: k = klayer(x,y,kref)
%
% Input: 
% [x,y]: meshgrid on (0,1)^2
%  kref: reference wavenumber (mid layer) 
%
% Output:
% k: matrix of layered wavenumbers on the meshgrid [x,y]
%
%  Author: Luis Garcia Ramos, 
%          Institut fur Mathematik, TU Berlin
%
%Version 0.1 - Nov 2015
%%%%%
%%
y1 = 0.2*x+0.2;
y2 = -0.2*x+0.8;

index1 = (y<=y1);
index2 = y>y1;
index3 = y>y2;

k = zeros(size(x));

k(index1) = 1.5*kref;
k(index2) = kref;
k(index3) = 2*kref;
end
