function [A] = helmholtz(k,np,bc)
%% HELMHOLTZ: Constructs matrices for the 1D Helmholtz problem.
%  Constructs the finite difference matrix corresponding
%  to the 1D Helmholtz problem discretization of -u''- k^2u = f
%  
%  When homogeneous Dirichlet boundary conditions are used, 
%  the boundary points are eliminated of the linear system.
%  In case of Sommerfeld boundary conditions the boundary points
%  are included.
%
%  Use: [A] = helmholtz(k,np,bc)
%
%  Input: 
%  k:      wavenumber
%  np:     number of interior discretization points
%  bc:     type of boundary conditions:      
%          'dir' for homogeneous dirichlet bc's
%          'som' for sommerfeld bc's 
%
%  Output:
%  A:      discretization matrix of Helmholtz problem
%
%  Author: Luis Garcia Ramos, 
%          Institut fur Mathematik, TU Berlin
%  Version 0.1 - Nov 2015
%
%  TO DO: Add 2D case
%         Include option for right hand side
%%%%%

%% Construction of 1D matrices
h  = 1/(np+1);            %gridsize

switch bc
    case 'dir'
        % Dirichlet 1D matrix (no boundary points)
        nv = np;                   %size of the linear system
        l  = ones(nv,1)*(-1/h^2);  %lower(=upper) diagonal
        d  = ones(nv,1)*(2/h^2); 
        A  = spdiags([l d l],[-1 0 1],nv,nv)- k^2*speye(nv); %Helmholtz matrix
        
        
    case 'som'
        %Sommerfeld 1D matrix (with boundary points)
        nv = np+2;
        d  = ones(np,1)*(2-k^2*h^2);    
        a  = (1-(k^2*h^2)/2+1i*k*h); % boundary conditions (second order)
        b  = (1-(k^2*h^2)/2-1i*k*h); 
        d  = [a; d; b];   
        u  = -ones(nv,1);      
        A  = 1/(h^2)*spdiags([u d u],[-1 0 1],np+2,np+2); %Helmholtz matrix
        
        %x  = 0:h:1;
        %b  = feval(f,x)'; b(1)=b(1)/2; b(length(b))=b(length(b))/2; %Right hand side
                        
        otherwise
            error('invalid boundary conditions')
end


end