function [A, sol,b] = helmholtz1d(f,k,np,bc,flag)
%% HELMHOLTZ1D: Direct solver for 1-D the Helmholtz equation.
%  Solves -u''-k^2*u=f with various boundary conditions
%  INPUT: 
%  f:      right-hand side (function handle)
%  k:      wavenumber of Helmholtz equation
%  np:     number of interior discretization points
%  bc:     type of boundary conditions:      
%          'dir' for homogeneous dirichlet bc's
%          'som' for sommerfeld bc's 
%  flag:   if flag==1 solve exactly and return solution
%
%  OUTPUT:
%  A:      discrete Helmholtz operator
%  b:      right hand side vector
%  sol:    solution of the linear system
%
%  AUTHOR: Luis Garcia Ramos, 
%        Institut fur Mathematik, TU Berlin
%  Version 0.1 - June 2015
%%%%%

%% Construction of 1D matrices
h  = 1/(np+1);            %gridsize

switch bc
    case 'dir'
        % Dirichlet 1D matrix (no boundary points)
        nv = np;                  %size of the linear system
        l  = ones(nv,1)*(-1/h^2); %lower(=upper) diagonal
        d  = ones(nv,1)*(2/h^2); 
        A  = spdiags([l d l],[-1 0 1],nv,nv)- k^2*speye(nv); %Helmholtz matrix
        size(A)
        x  = h:h:1-h;
        b  = feval(f,x)'; %right hand side
        
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
        b  = feval(f,x)'; b(1)=b(1)/2; b(length(b))=b(length(b))/2; %Right hand side
                        
        otherwise
            error('invalid boundary conditions')
end

% solve directly the linear system
sol = zeros(size(b));
if flag == 1
    sol = A\b;  % Solution 
end

end
