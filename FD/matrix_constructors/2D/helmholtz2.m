function [A] = helmholtz2(k,eps,npx,npy,bc)
%  HELMHOLTZ2: Constructs matrices for the 2D Helmholtz
%  and shifted Laplace problems.
%  Constructs the finite difference matrix corresponding
%  to the discretization of div(grad u)- (k^2+i*eps)u = f in (0,1)x(0,1)
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
%Version 1.0 - Jul 2016
%%%%%

%% Construction of 2D matrix

hx  = 1/(npx+1);  %gridsize in x-direction
hy  = 1/(npy+1);  %gridsize in y-direction

switch bc
    case 'dir'        
        l    = ones(npx,1)*(-1/hx^2); 
        d    = ones(npx,1)*(2/hx^2); 
        Dx   = spdiags([l d l],[-1 0 1],npx,npx); %1D Neg Laplacian, x direction
        l    = ones(npy,1)*(-1/hy^2);
        d    = ones(npy,1)*(2/hy^2);
        Dy   = spdiags([l d l],[-1 0 1],npy,npy); %1D Neg Laplacian, y direction
        NLap = kron(Dx,speye(npy))+kron(speye(npx),Dy); %2D Negative Laplacian
        A    = NLap - (k^2+1i*eps)*speye(length(NLap));
        
    case 'som'
        %2D matrix with Sommerfeld bc's (with boundary points)
        nx = npx+2;
        ny = npy+2;
        
        %Interior points (i*hx,j*hy)
        l    = ones(nx,1)*(-1/hx^2); 
        d    = ones(nx,1)*(2/hx^2); 
        Dx   = spdiags([l d l],[-1 0 1],nx,nx); %1D Laplacian, x direction
        l    = ones(ny,1)*(-1/hy^2);
        d    = ones(ny,1)*(2/hy^2);
        Dy   = spdiags([l d l],[-1 0 1],ny,ny); %1D Laplacian, y direction
        NLap = kron(Dx,speye(ny))+kron(speye(nx),Dy); %2D Negative Laplacian
        A    = NLap - (k^2+1i*eps)*speye(nx*ny);
        
        %Noncorner boundary points
        
        %West boundary: (0,j*hy)
         for j = 1:npy
              n = j*(npx+2)+1; %changing from (0,j) to global index
              A(n,:)       = zeros(nv,1)'; %clear row
              A(n,n)       = (-k^2-1i*eps-2*1i*k/hx+2/hx^2+2/hy^2);
              A(n,n+1)     = -2/hx^2;
              A(n,n+npx+2) = -1/hy^2;
              A(n,n-npx-2) = -1/hy^2;   
         end
         
             
         %South boundary: (ihx,0)
         for i=1:npx
              %k=k;
              n=i+1; %changing from (i,0) to global index
              A(n,:)       = zeros(nv,1)';
              A(n,n)       = -k^2-1i*eps-2*1i*k/hy+2/hx^2+2/hy^2;
              A(n,n+1)     = -1/hx^2;
              A(n,n-1)     =  -1/hx^2;
              A(n,n+npx+2) = -1/hy^2;  
         end
         
          
         %East boundary: (1,j*hy)
         for j=1:npy
           % k = k;
             n = (npx+2)*(j+1); %changing from (npx+1,j) to global index
             A(n,:)       = zeros(nv,1)';
             A(n,n)       = (-k^2-1i*eps-2*1i*k/hx+2/hx^2+2/hy^2);
             A(n,n-1)     = -1/hx^2;
             A(n,n-npx-2) = -1/hy^2;  
             A(n,n+npx+2) = -1/hy^2;
         end

         %North boundary: (i*hx,1)
         for i=1:npx
             n = (npx+2)*(npy+1)+(i+1);
             A(n,:)       = zeros(nv,1)';
             A(n,n)       = (-k^2-1i*eps-2*1i*k/hy+2/hx^2+2/hy^2);
             A(n,n-1)     =  -1/hx^2;
             A(n,n+1)     =  -1/hx^2;  
             A(n,n-npx-2) = -2/hy^2;
         end
         
          %SW Corner: (0,0)
          %k = k;
          A(1,:)       = zeros(nv,1)';
          A(1,1)     = (-k^2-1i*eps-2*1i*k/hy-2*1i*k/hx+2/hx^2+2/hy^2);
          A(1,2)     = -2/hx^2;
          A(1,npx+3) = -2/hy^2;
          
          %SE Corner: (1,0)
          %k = k;
          n = npx+2;
          A(n,:)    = zeros(nv,1)';
          A(n,n)    = -k^2-1i*eps-2*1i*k/hy-2*1i*k/hx+2/hx^2+2/hy^2;
          A(n,n-1)  = -2/hx^2;
          A(n,2*n)  = -2/hy^2;
          
          %NW Corner: (0,1)
          %k = k;
          n = (npx+2)*(npy+1)+1;
          A(n,:)         = zeros(nv,1)';
          A(n,n)         = -k^2-1i*eps-2*1i*k/hy-2*1i*k/hx+2/hx^2+2/hy^2;
          A(n,n+1)       = -2/hx^2;
          A(n,n-(npx+2)) = -2/hy^2;

          
          %NE Corner: (1,1)
          %k = k;
          n = (npx+2)*(npy+2);
          A(n,:)         = zeros(nv,1)';
          A(n,n)         = -k^2-1i*eps-2*1i*k/hy-2*1i*k/hx+2/hx^2+2/hy^2;
          A(n,n-1)       = -2/hx^2;
          A(n,n-(npx+2)) = -2/hy^2;

        otherwise
            error('invalid boundary conditions')
end

end



