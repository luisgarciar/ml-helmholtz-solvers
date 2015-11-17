function [S, sol ] = SLaplace_1D(f,sigma,beta1,beta2,flag)
%SLAPLACE_1D: Direct solver for shifted Laplacian system with parameters 
%              sigma, beta1, beta2
%
%INPUT: 
%  f:      right-hand side function without values on the boundary 
%  sigma:  wavenumber of Helmholtz equation
%  beta1:  real shift 
%  beta2:  complex shift


%OUTPUT:
%  S:      shifted Laplace operator without boundary conditions
%  sol:    exact solution without boundary conditions

nl   = length(f)+1;     %number of interval 
nv   = length(f);
h    = 1/nl;         % gridsize
beta = complex(beta1,-beta2) % shift parameter

d   = ones(nv,1)*(-sigma*(beta)+(2/h^2));
l   = ones(nv,1)*(-1/h^2); 
S   = spdiags([l d l],[-1 0 1],nv,nv); % Discrete operator of interior gridpoints
sol = zeros(size(f));

if flag==1
    sol = S\f;  % Solution without boundary conditions (zeros from Dirichlet)
end

end

