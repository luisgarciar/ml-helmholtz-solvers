function [A,sol,b] = helmholtzBC(f,k,np,bc,dim,flag)
%% HELMHOLTZBC: Direct solver for the Helmholtz equation.
%Solves -u''+k^2*u=f with various boundary conditions

%INPUT: 
%  f:      right-hand side function handle 
%  k:      wavenumber of Helmholtz equation
%  nd:     number of interior discretization points in 1D           
%  bc:     type of boundary conditions:      
%          'dir' for homogeneous dirichlet bc's
%          'som' for sommerfeld bc's 
%  dim:    problem dimension (1D or 2D)
%  flag:   if flag==1 solve exactly and return solution

%OUTPUT:
%  A:      Discrete Helmholtz operator 
%  sol:    Solution of the linear system
%
% Note: 
% For homogeneous Dirichlet problems, the linear system is obtained
% after elimination of the boundary conditions.
% For Sommerfeld problems, the boundary points are not eliminated.



%% Construction of 1D matrices
n  = np + 1; 
h  = 1/n;         %gridsize
nv = np;       %size of the linear system

% Dirichlet 1D matrix (no boundary points)
l  = ones(nv,1)*(-1/h^2); %lower(=upper) diagonal
d  = ones(np,1)*(2/h^2-k^2);
Ad_1 = spdiags([l d l],[-1 0 1],nv,nv);


%Sommerfeld 1D matrix with boundary points
%See Elman, O'Leary, Numer. Math. Vol. 83, Issue 2, p. 231-257, 1999)
d     = ones(np,1)*(2/h^2-k^2); 
gamma = (1i*k-1/h); 
d     =  [gamma; d; gamma];  %diagonal
u     =  ones(n,1)*(-1/h^2); 
u     =  [1/h ; u];          %upp diag
l     =  [u   ; 1/h];        %lower diag  
As_1  = spdiags([l d u],[-1 0 1],nv,nv); 


%% Construction of 2D matrices
%(See Elman, O'Leary, Numer. Math., 1999)

% Homogeneous Dirichlet bc's
d    = ones(np,1)*(2/h^2); 
l    = ones(np,1)*(-1/h^2); %upper diagonal
B    = spdiags([l d l],[-1 0 1],np,np);
Ad_2 = kron(B, speye(np)) + kron(speye(np), B) - k^2*speye(np^2) ;

% Sommerfeld bc's
A0     = As_1 + k^2*speye(np);
As_2   = kron(speye(np),A0)   + kron(As_1,speye(np));

%% Solution of the linear system
% Construction of mesh and right hand side b (point source)for 1D
 if dim==1
 x = h:h:1-h;
 b = feval(f,x);
 b = [0;b;0];
 end

% Construction of mesh and right hand side b (point source) for 2D
 if dim==2
[x,y] = meshgrid(h:h:1-h);
 b    = feval(f,x,y);
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
    sol = A\b;  % Solution without boundary conditions    
end

end