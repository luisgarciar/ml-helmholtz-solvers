function [A] = helmholtz2_ord1(k,eps,npx,npy,bc)
%% HELMHOLTZ2_ORD1: Constructs matrices for the 2D Helmholtz and shifted Laplace problems.
%  Constructs the finite difference matrix corresponding
%  to the discretization of 
%  div(grad u)- (k^2+i*eps)u = f in (0,1)x(0,1)
%
% With boundary conditions
%  u  =   0 on boundary  (Dirichlet)
%  or
%  du/dn-iku = 0 (Sommerfeld)
%
% The boundary condition is discretized with 
% one-sided differences, and the discretization
% is O(h).
%
% Reference:
% Ernst, O. and Golub, G., A domain decomposition approach to solving 
% the Helmholtz equation with a radiation boundary condition (1994)
%
%  Use: [A] = helmholtz2(k,eps,npx,npy,bc)
%
%  When homogeneous Dirichlet boundary conditions are used, 
%  the boundary points are eliminated of the linear system.
%  In case of Sommerfeld boundary conditions the boundary points
%  are included.
%
%  Note: The imaginary shift eps is not multiplied by k^2 as in the papers
%        of Erlangga et al.
%
%  For convergence of multigrid applied to the shifted Laplacian
%  one needs eps~k^2 (typically eps=0.5*k^2)
%
%  For a number of iterations of preconditioned GMRES independent of k
%  one needs eps~k (but multigrid will fail for the shifted Laplacian)
%
%  Input: 
%  k:      real wavenumber of Helmholtz equation 
%  eps:    imaginary part of the shift (eps=0 for pure Helmholtz problem)
%  npx:    number of interior discretization points in the x-direction
%  npy:    number of interior discretization points in the y-direction
%  bc:     type of boundary conditions:      
%          'dir' for homogeneous dirichlet bc's
%          'som' for 1st order sommerfeld bc's 
%
%  Output:
%  A:      discretization matrix of Helmholtz problem
%          size(A) = (npx,npy)     for Dirichlet  problems    
%                  = (npx+2,npy+2) for Sommerfeld problems
%
%  Author: Luis Garcia Ramos, 
%          Institut fur Mathematik, TU Berlin
%
%Version 2.0 - Jul 2016
%%%%%

%% Construction of 2D matrix

np  = max(npx,npy);
h   = 1/(np+1); 
hx  = 1/(npx+1);  %gridsize in x-direction
hy  = 1/(npy+1);  %gridsize in y-direction

switch bc
    case 'dir'                
          np = npx*npy;
          W = -ones(np,1)/hx^2;  E=W; %Dxx
          N = -ones(np,1)/hy^2; S=N; %Dyy
          C = 2*ones(np,1)/hx^2+2*ones(np,1)/hy^2;
          A = spdiags([S W C E N],[-npx -1 0 1 npx],np, np)-(k^2+1i*eps)*speye(np);
          
          for i=1:(npy-1)      %Modify points closest to east and west boundaries
              ii=npx*(i-1)+npx;
              A(ii,ii+1) = 0;
              A(ii+1,ii) = 0;
          end
                    
    case 'som'        
        %2D matrix with Sommerfeld bc's (with boundary points)
        npts = max(npx,npy);
        h    = 1/(npts+1); 
        np   = (npts+2)^2;
        
        %For the construction see the paper in the references
        d = 4-(k^2+1i*eps)*h^2-1i*k*h; 
        T = gallery('tridiag',npts+2,-1,d,-1);
        T(1,1) = 3-(k^2+1i*eps)*h^2-1i*k*h;
        T(npts+2,npts+2)= 3-(k^2+1i*eps)*h^2-1i*k*h;
        alpha = -1-1i*k*h;          
       
        J = gallery('tridiag',npts+2,1,0,1);
        V = sparse(npts+2,npts+2);
        V(1,1) = 1; V(npts+2,npts+2) = 1;
        I = speye(npts+2);
        
        A = kron(I,T)+kron(alpha*V-J,I);
        
    otherwise
        error('invalid boundary conditions')
end

end



