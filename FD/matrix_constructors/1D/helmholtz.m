function [A] = helmholtz(k,eps,np,bc)
%% HELMHOLTZ: Constructs matrices for the 1D Helmholtz problem.
%  Constructs the finite difference matrix corresponding
%  to the discretization of the 1D Helmholtz/shifted Laplace problem
%
%       -u''- k^2u = f in (0,1)
%        u(0)=0, u(1)=0, or
%        u'(0)-iku(0)=0, u'(1)+iku(1)=1
%  
%  Use: [A] = helmholtz(k,eps,np,bc)

% Note (For shifted Laplacian problems): 
% The imaginary shift is not multiplied by the term k^2
% as in the original publications of Erlangga et. al.
% For a convergent multigrid solver one needs eps~k^2 (eps=0.5k^2)
% For wave-independent GMRES convergence one needs eps~k 
% (but the MG solver will not converge)
%
%
%  Input: 
%  k:      wavenumber
%  eps:    imaginary part of shifted Laplacian
%  np:     number of interior discretization points
%  bc:     type of boundary conditions:      
%          'dir' for homogeneous dirichlet bc's
%          'som' for sommerfeld bc's 
%
%  Output:
%  A:      discretization matrix of Helmholtz/shifted Laplace problem
%
%  Author: Luis Garcia Ramos, 
%          Institut fur Mathematik, TU Berlin
%  Version 1.0 - Jul 2016
%
%%%%%

%% Construction of 1D matrices
h  = 1/(np+1);            %gridsize

switch bc
    case 'dir'
        % Dirichlet 1D matrix (no boundary points)
        nv = np;                   %size of the linear system
        l  = ones(nv,1)*(-1/h^2);  %lower(=upper) diagonal
        d  = ones(nv,1)*(2/h^2); 
        A  = spdiags([l d l],[-1 0 1],nv,nv)- (k^2+1i*eps)*speye(nv); %Helmholtz matrix
        
    case 'som'
        %Sommerfeld 1D matrix (with boundary points)
        nv = np+2;
        d  = ones(np,1)*(2-(k^2+1i*eps)*h^2);    
        a  = (1-(k^2+1i*eps)*h^2/2+1i*k*h); % boundary conditions (second order)
        b  = (1-(k^2+1i*eps)*h^2/2-1i*k*h); 
        d  = [a; d; b];   
        u  = -ones(nv,1);     
        A  = 1/(h^2)*spdiags([u d u],[-1 0 1],np+2,np+2); %Helmholtz matrix
        
                  
        otherwise
            error('invalid boundary conditions')
end


end