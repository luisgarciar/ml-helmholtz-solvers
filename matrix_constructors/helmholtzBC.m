function [A,sol,b] = helmholtzbc(dim,bc,np,k,f,flag)
%% HELMHOLTZBC: Direct solver for the Helmholtz equation.
%  Solves -div(grad u)+k^2*u=f with various boundary conditions
%
%  USAGE: [A,sol,b] = helmholtzbc(dim,bc,np,k,f,flag)
% 
%  INPUT: 
%  dim:    problem dimension (1D or 2D)
%  bc:     type of boundary conditions:     
%          'dir' for homogeneous dirichlet bc's
%          'som' for first order sommerfeld bc's
%  np:     number of interior discretization points in 1D
%  k:      wavenumber of Helmholtz equation (only constant case)           
%  f:      right-hand side (function handle) 
%  flag:   if flag==1 solve exactly and return solution
%
% OUTPUT:
%  A:      discrete Helmholtz operator
%  b:      right hand side vector
%  sol:    solution of the linear system
%
% Note: 
% For homogeneous Dirichlet problems, the linear system is obtained
% after elimination of the boundary conditions.
% For Sommerfeld problems, the boundary points are not eliminated.
%
% Author: Luis Garcia Ramos, 
%         Institut fur Mathematik, TU Berlin
% Version 0.2 - Nov 2015
%
% TO DO: Fix 2D case- Sommerfeld BC's
%        Add non-constant wavenumber

%% Checking grid size and wavenumber
n  = np + 1; 
h  = 1/n;      %gridsize 

if np < ceil(5*k/pi)-1
    %disp('The number of points per wavelength is too small. Setting to minimum = 10ppw')
    %np = round(5*k/pi)-1;
end

%% Construction of the matrix and right hand side
switch dim
    case 1
    switch bc
         case 'dir'
             % Dirichlet 1D matrix (no boundary points)
             nv = np;   %size of the linear system
             l  = ones(nv,1)*(-1/h^2); %lower(=upper) diagonal
             d  = ones(nv,1)*(2/h^2); 
             A  = spdiags([l d l],[-1 0 1],nv,nv)- k^2*speye(nv); %Helmholtz matrix             
             
             x  = h:h:1-h;
             b  = feval(f,x,h); %right hand side
 

        case 'som'
            %Sommerfeld 1D matrix (with boundary points)
            nv = np+2;
            d  = ones(np,1)*(2-k^2*h^2);    
            a  = (1-(k^2*h^2)/2+1i*k*h); % boundary conditions
            b  = (1-(k^2*h^2)/2-1i*k*h); 
            d  = [a; d; b];   
            u  = -ones(nv,1);      
            A  = 1/(h^2)*spdiags([u d u],[-1 0 1],np+2,np+2); %Helmholtz matrix
            
            x  = 0:h:1;
            b  = feval(f,x,h); b(1)=b(1)/2; b(length(b))=b(length(b))/2; %Right hand side
                        
        otherwise
            error('invalid boundary conditions')
    end
    
    case dim==2
        switch bc
            case 'dir'
                % Dirichlet 2D matrix (no boundary points)
                d    = ones(np,1)*(2/h^2); 
                l    = ones(np,1)*(-1/h^2); %upper diagonal
                B    = spdiags([l d l],[-1 0 1],np,np);
                A    = kron(B, speye(np)) + kron(speye(np), B) - k^2*speye(np^2);
                
                [x,y] = meshgrid(h:h:1-h);
                   b  = feval(f,x,y);
                   b  = reshape(b',[np^2,1]);
                
            case 'som'
                %Sommerfeld 2D matrix (with boundary points) (TO BE FIXED)
                nv    = np+2;
                d     = ones(np,1)*(2-k^2*h^2);    
                a     = (1-(k^2*h^2)/2+1i*k*h); %boundary conditions
                b     = (1-(k^2*h^2)/2-1i*k*h); 
                d     = [a; d; b];   
                u     = -ones(nv,1);      
                A  = 1/(h^2)*spdiags([u d u],[-1 0 1],np+2,np+2); %Helmholtz matrix
                
                % Construction of right hand side (TO BE FIXED)
                [x,y] = meshgrid(0:h:1);
                   b  = feval(f,x,y);
                   b  = reshape(b',[nv,1]);
                   
            otherwise
                error('invalid boundary conditions')
        end
        
    otherwise
        error('invalid dimension')
end

%% Solution of the linear system
A   = sparse(A);
% solve directly the linear system
sol = zeros(size(b));
if flag == 1
    sol = A\b;  % Solution 
end

end