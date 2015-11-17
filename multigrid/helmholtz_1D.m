function [A, sol ] = helmholtz_1D(f,sigma,flag)
%HELMHOLTZ_1D: Direct solver for the 1-D Helmholtz equation.
%Solves -u''+sigma*u=f with Dirichlet homogeneous boundary conditions
%
%INPUT: 
%  f:      right-hand side function without values on the boundary 
%  sigma:  wavenumber of Helmholtz equation
%  flag:   if flag==1 solve exactly and return solution


%OUTPUT:
%  A:      Discrete Helmholtz operator without boundary conditions
%  sol:    exact solution without boundary conditions

nl=length(f)+1;     %number of interval 
nv=length(f)
h   = 1/nl;         % gridsize

d   = ones(nv,1)*(-sigma+2/h^2);
l   = ones(nv,1)*(-1/h^2); 
A   = spdiags([l d l],[-1 0 1],nv,nv); % Discrete operator of interior gridpoints
size(A)
size(f)

sol = zeros(size(f));

if flag==1
    sol = A\f;  % Solution without boundary conditions (zeros from Dirichlet)
end

end

