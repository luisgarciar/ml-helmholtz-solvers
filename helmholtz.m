function [A, sol,b] = helmholtz(f,k,nd,bc,dim,flag)
%% HELMHOLTZ: Direct solver for the Helmholtz equation.
%Solves -u''+k^2*u=f with various boundary conditions

%INPUT: 
%  f:      right-hand side function
%  k:      wavenumber of Helmholtz equation
%  nd:     number of interior discretization points
%  bc:     type of boundary conditions:      
%          'dir' for homogeneous dirichlet bc's
%          'som' for sommerfeld bc's 
%  dim:    problem dimension (1D or 2D)
%  flag:   if flag==1 solve exactly and return solution


%OUTPUT:
%  A:      Discrete Helmholtz operator with homogeneous Dirichlet boundary conditions after elimination of boundary conditions
%  sol:    Solution without boundary conditions

%% Construction of 1D matrices
h  = 1/(nd+1);         %gridsize
l  = ones(nd,1)*(-1/h^2); %upper diagonal
d  = ones(nd,1)*(2/h^2-k^2); 

% Dirichlet 1D matrix (only interior points)
d      = ones(nd,1)*(2/h^2); 
Ad_1   = spdiags([l d l],[-1 0 1],nd,nd)- k^2*speye(nd); 

% Sommerfeld 1D matrix (after elimination of bc's, see Elman, O'Leary, Numer. Math. 1999)
d      = ones(nd,1)*(2/h^2-k^2); 
gamma  = 2/h^2-k^2-(1+1i*k*h)/(h^2*(1+k^2*h^2));
d(1)   = gamma;
d(nd)  = gamma;

As_1   = spdiags([l d l],[-1 0 1],nd,nd); 
%% Construction of 2D matrices

% Homogeneous Dirichlet bc's
d    = ones(nd,1)*(2/h^2); 
l    = ones(nd,1)*(-1/h^2); %upper diagonal
B    = spdiags([l d l],[-1 0 1],nd,nd);
Ad_2 = kron(B, speye(nd)) + kron(speye(nd), B) - k^2*speye(nd^2) ;

% Sommerfeld bc's
A0     = As_1 + k^2*speye(nd);
As_2   = kron(speye(nd),A0)   + kron(As_1,speye(nd));

%% Solution of the linear system
% Construction of mesh and right hand side b for 1D
 if dim==1
 x = h:h:1-h;
 b = feval(f,x);
 end

% Construction of mesh and right hand side b for 2D
 if dim==2
[x,y] = meshgrid(h:h:1-h);
 b    = feval(f,x,y);
 b    = reshape(b',[nd^2,1]);
 end

% Set the matrix of the linear system according to bc's and dim
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
    sol = A\b;  % Solution without boundary conditions    
end

end