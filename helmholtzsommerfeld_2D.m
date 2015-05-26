function [A, sol] = helmholtzsommerfeld_2D(f,k,nd,flag)
% HELMHOLTZSOMMERFELD_2D: Direct solver for the 2-D Helmholtz equation.
% Solves -u''+ ku = f with Sommerfeld boundary conditions
% on a uniform discretization of the square [0,1]x[0,1] 

% INPUT: 
%  f:      right-hand side function without values on the boundary
%  k:      wavenumber of Helmholtz equation
%  nd:     number of interior discretization points in one direction
%  flag:   if flag == 1 solve exactly and return solution

% OUTPUT:
%  A:      Discrete Helmholtz operator after elimination of boundary conditions
%  sol:    Solution without boundary conditions

if (isequal(nd^2, length(f)) ~= 1)
    error('Invalid right-hand side: size must correspond to a square grid')
end

%% Construction of the 1D Sommerfeld matrix
nv = length(f);        %number of variables
h  = 1/(nd+1);         % gridsize

% Diagonal of the matrix (after elimination of bc's)
d       = ones(nv,1)*(2/h^2-k^2); 
gamma   = 2/h^2-k^2-(1+1i*k*h)/(h^2*(1+k^2*h^2));
d(1)    = gamma;
d(nv)   = gamma;

l       = ones(nv,1)*(-1/h^2); %upper diagonal

% 1-D matrix
A   = spdiags([l d l],[-1 0 1],nv,nv); 

%% Construction of the 2D Sommerfeld matrix (See Elman, O'Leary, Numer. Math. 1999)
A0  = A + k^2*eye(nv^2);
A   = kron(A0,eye(nv^2)) + kron(eye(nv^2),A);
%A   = sparse(A);

%% Solution of the linear system,
sol = zeros(nv);
if flag == 1
    sol = A\f;  % Solution without boundary conditions    
end



