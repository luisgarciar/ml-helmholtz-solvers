function [A, sol,b] = helmholtz1(f,k,np,bc,dim,flag)
%% HELMHOLTZ1: Direct solver for the Helmholtz equation.
%  Solves -u''-k^2*u=f with various boundary conditions
%  INPUT: 
%  f:      right-hand side (function handle)
%  k:      wavenumber of Helmholtz equation
%  np:     number of interior discretization points
%  bc:     type of boundary conditions:      
%          'dir' for homogeneous dirichlet bc's
%          'som' for sommerfeld bc's 
%  dim:    problem dimension (1D or 2D)
%  flag:   if flag==1 solve exactly and return solution
%
% OUTPUT:
%  A:      Discrete Helmholtz operator with homogeneous Dirichlet boundary conditions after elimination of boundary conditions
%  sol:    Solution without boundary conditions
%
% AUTHOR: Luis Garcia Ramos, 
%        Institut fur Mathematik, TU Berlin
% Version 0.1 - June 2015

%% Construction of 1D matrices
h  = 1/(np+1);            %gridsize
l  = ones(np,1)*(-1/h^2); %upper diagonal

% Homogeneous Dirichlet boundary conditions
nv      = np;
d       = -2*ones(nv,1); 
l       = ones(nv,1); %upper diagonal

Lapl1d  = spdiags([l d l],[-1 0 1],nv,nv); %1D Discrete Laplacian 
Ad_1   = -(1/h^2)*Lapl1d - k^2*speye(nv);   %1D Discrete Helmholtz Operator

%CHECK SOMMERFELD MATRIX!
%Sommerfeld boundary conditions
%(See Elman, O'Leary, Numer. Math. Vol. 83, Issue 2, p. 231-257, 1999)
d      = ones(np,1)*(2/h^2-k^2); 
gamma  = 2/h^2-k^2-(1+1i*k*h)/(h^2*(1+k^2*h^2));
d(1)   = gamma;
d(np)  = gamma;

As_1   = spdiags([l d l],[-1 0 1],np,np); 


%% Construction of 2D matrices
%(See Elman, O'Leary, Numer. Math., 1999)
nv = np^2;

% Homogeneous Dirichlet boundary conditions
Lapl2d =  kron(Lapl1d, speye(np)) + kron(speye(np), Lapl1d);
Lapl2d = (1/h^2)*Lapl2d; 

Ad_2   =  -Lapl2d - k^2*h^2*speye(nv);  %2D Discrete Helmholtz Operator


%CHECK SOMMERFELD MATRIX!
% Sommerfeld bc's
A0     = As_1 + k^2*speye(np);
As_2   = kron(speye(np),A0)   + kron(As_1,speye(np));

%% Solution of the linear system
% Construction of mesh and right hand side b (point source)for 1D
 if dim==1
 x = h:h:1-h;
 b = feval(f,x);
 end

% Construction of mesh and right hand side b (point source) for 2D
 if dim==2
 [x,y] = meshgrid(h:h:1-h);
 b    = h^2*feval(f,x,y);
 b    = reshape(b',[np^2,1]);
 end

% Set the matrix of the linear system according to boundary conditions
% and dimension of the problem
switch bc
    case 'dir'
        if dim==1
            A = Ad_1;
        elseif dim==2
            A = Ad_2;
        else
            error('invalid dimension') 
        end
        
    case 'som'
        if dim==1
            A = As_1;
        elseif dim==2
            A = As_2;
        else
            error('invalid dimension'); 
        end
        
    otherwise
        error('invalid boundary conditions'); 
end

A   = sparse(A);
% solve directly the linear system
sol = zeros(size(b));
if flag == 1
    sol = A\b;  % Solution 
end

end