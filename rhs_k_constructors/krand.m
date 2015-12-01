function k = krand(x,y,kmin,kmax)
%% KRAND Constructs a matrix of random wavenumbers for the Helmholtz problem
%        
% Use: f = krand(x,y,kmin,kmax)
%
% Input: 
% [x,y]: meshgrid on (0,1)^2
%  kmin, kmax: minimum (resp. maximum) wavenumber
%
% Output:
% k: matrix of (uniformly distributed) random wavenumbers on the meshgrid [x,y]
%
%  Author: Luis Garcia Ramos, 
%          Institut fur Mathematik, TU Berlin
%%
%Version 0.1 - Nov 2015
%%%%%

%%
assert(isequal(size(x),size(y)),'Inconsistent input size'); 
k = kmin+(kmax-kmin).*rand(size(x));

end
