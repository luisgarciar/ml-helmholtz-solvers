function [A, sol ] = helmholtzsommerfeld_1D(f,k,flag)
%% HELMHOLTZSOMMERFELD_1D: Direct solver for the 1-D Helmholtz equation.
%Solves -u''+sigma*u=f with Sommerfeld boundary conditions
%
%INPUT: 
%  f:      right-hand side function 
%  k:      wavenumber of Helmholtz equation
%  flag:   if flag==1 solve exactly and return solution

%OUTPUT:
%  A:      Discrete Helmholtz operator after elimination of boundary conditions
%  sol:    Solution without boundary conditions

%% 
nl = length(f)+1;     %number of interval 
nv = length(f)
h  = 1/nl;         % gridsize

%diagonal of the matrix (after elimination of bc's)
d       = ones(nv,1)*(2/h^2-k^2);
gamma   = 2/h^2-k^2-(1+1i*k*h)/(h^2*(1+k^2*h^2))
d(1)    = gamma
d(nv)   = gamma;
l       = ones(nv,1)*(-1/h^2); 

A   = spdiags([l d l],[-1 0 1],nv,nv); % Discrete operator of interior gridpoints

sol = zeros(size(f));

if flag==1
    sol = A\f;  % Solution without boundary conditions (zeros from Dirichlet)    
end