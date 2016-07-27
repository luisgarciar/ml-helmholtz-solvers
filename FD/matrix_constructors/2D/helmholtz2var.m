function [A] = helmholtz2var(kvar,epsvar,npx,npy,bc)
%% HELMHOLTZ2VAR: Constructs matrices for the 2D Helmholtz 
%  and shifted Laplace problems with non-constant wavenumbers.
%
%  Constructs the finite difference matrix corresponding
%  to the discretization of 
%
%  div(grad u)-(kvar^2 + i*epsvar)*u = f in (0,1)x(0,1)
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
%  When homogeneous Dirichlet boundary conditions are used, 
%  the boundary points are eliminated of the linear system.
%  In case of Sommerfeld boundary conditions the boundary points
%  are included.
%
%  The discretization of Sommerfeld BC's is of second order
%
%  Use: [A] = helmholtz2(kvar,epsvar,npx,npy,bc)
%
%  Input: 
%  kvar:   wavenumber (function handle, k depends on (x,y))
%  eps:    imaginary shift (function handle, k depends on (x,y))
%  npx:    number of interior discretization points in the x-direction
%  npy:    number of interior discretization points in the y-direction
%  bc:     type of boundary conditions:      
%          'dir' for homogeneous dirichlet bc's
%          'som' for 1st order sommerfeld bc's 
%
%  Output:
%  A:      discretization matrix of Helmholtz problem
%          size(A)=(npx*npy,npx*npy) for Dirichlet problems
%                 =((npx+2)*(npy+2),(npx+2)*(npy+2)) for Sommerfeld prob.
%          
%  Author: Luis Garcia Ramos, 
%          Institut fur Mathematik, TU Berlin
%
%Version 1.0 - Nov 2016
%%%%%
%% Construction of 2D matrix
hx  = 1/(npx+1);            %gridsize
hy  = 1/(npy+1);            %gridsize

switch bc
    case 'dir'
        % Dirichlet 2D matrix (no boundary points)        
        nv = npx*npy;
        v  = ones(nv,1);
        N  = (-1/hy^2)*v; S=N;
        W  = (-1/hx^2)*v; E=W;
        C  = (2/hx^2+2/hy^2)*v;
        
        %Create vector of wavenumbers
        [x,y] = meshgrid(hx:hx:1-hx,hy:hy:1-hy);
         kk   = feval(kvar,x,y);
         kk   = reshape(kk',[nv,1]);
         eps  = feval(epsvar,x,y);
         eps  = reshape(eps',[nv,1]);
        
        %A: 2D Helmholtz matrix
        A = spdiags([S W C E N],[-npx -1 0 1 npx], nv, nv);
        
        for i=1:(npy-1)       %Modify points closest to east and west boundaries
              ii=npx*(i-1)+npx;
              A(ii,ii+1) = 0;
              A(ii+1,ii) = 0;
        end
            
        A = A-(spdiags(kk,0,nv,nv).*spdiags(kk,0,nv,nv))...
              -1i*spdiags(eps,0,nv,nv);
         
    case 'som'        
        %2D matrix with Sommerfeld bc's (with boundary points)
        nv = (npx+2)*(npy+2);
        
        %Interior points (i*hx,j*hy)
        v = ones(nv,1);
        N = -1/hy^2*v; S=N;
        W = -1/hx^2*v; E=W;
        C = (2/hy^2+2/hy^2)*v;
        
        %A: 2D Helmholtz matrix (without the k^2 term)
        A = spdiags([S W C E N],[-(npx+2) -1 0 1 (npx+2)], nv, nv);         
        
        %Noncorner boundary points    
        %West boundary: (0,j*hy)
         for j = 1:npy
              k = feval(kvar,0,j*hy);
              n = j*(npx+2)+1; %changing from (0,j) to global index
              A(n,:)       = zeros(nv,1)'; %clear row
              A(n,n)       = (-2*1i*k/hx+2/hx^2+2/hy^2);
              A(n,n+1)     = -2/hx^2;
              A(n,n+npx+2) = -1/hy^2;
              A(n,n-npx-2) = -1/hy^2;   
         end
         
         %South boundary: (ihx,0)
         for i=1:npx
              k = feval(kvar,i*hx,0);
              n=i+1; %changing from (i,0) to global index
              A(n,:)       = zeros(nv,1)';
              A(n,n)       = (-2*1i*k/hy+2/hx^2+2/hy^2);
              A(n,n+1)     = -1/hx^2;
              A(n,n-1)    =  -1/hx^2;
              A(n,n+npx+2) = -1/hy^2;  
         end
           
         %East boundary: (1,j*hy)
         for j=1:npy
             k = feval(kvar,1,j*hy);
             n = (npx+2)*(j+1); %changing from (npx+1,j) to global index
             A(n,:)       = zeros(nv,1)';
             A(n,n)       = (-2*1i*k/hx+2/hx^2+2/hy^2);
             A(n,n-1)     = -1/hx^2;
             A(n,n-npx-2) = -1/hy^2;  
             A(n,n+npx+2) = -1/hy^2;
         end

         %North boundary: (i*hx,1)
         for i=1:npx
             k = feval(kvar,i*hx,1);
             n = (npx+2)*(npy+1)+(i+1);
             A(n,:)         = zeros(nv,1)';
             A(n,n)         = -2*1i*k/hy+2/hx^2+2/hy^2;
             A(n,n-1)       =  -1/hx^2;
             A(n,n+1)       =  -1/hx^2;  
             A(n,n-npx-2)   = -2/hy^2;
         end
         
          %SW Corner: (0,0)
          k = feval(kvar,0,0);
          A(1,:)      = zeros(nv,1)';
          A(1,1)     = -2*1i*k/hy-2*1i*k/hx+2/hx^2+2/hy^2;
          A(1,2)     = -2/hx^2;
          A(1,npx+3) = -2/hy^2;
          
          %SE Corner: (1,0)
          n = npx+2;
          k = feval(kvar,1,0);
          A(n,:)    = zeros(nv,1)';
          A(n,n)    = -2*1i*k/hy-2*1i*k/hx+2/hx^2+2/hy^2;
          A(n,n-1)  = -2/hx^2;
          A(n,2*n)  = -2/hy^2;
          
          %NW Corner: (0,1)
          k = feval(kvar,0,1);
          n = (npx+2)*(npy+1)+1;
          A(n,:)         = zeros(nv,1)';
          A(n,n)         = -2*1i*k/hy-2*1i*k/hx+2/hx^2+2/hy^2;
          A(n,n+1)       = -2/hx^2;
          A(n,n-(npx+2)) = -2/hy^2;

          
          %NE Corner: (1,1)
          k = feval(kvar,1,1);
          n = (npx+2)*(npy+2);
          A(n,:)         = zeros(nv,1)';
          A(n,n)         = -2*1i*k/hy-2*1i*k/hx+2/hx^2+2/hy^2;
          A(n,n-1)       = -2/hx^2;
          A(n,n-(npx+2)) = -2/hy^2;
          
          %Adding the effect of the wavenumber k and the imaginary shift on A
          %Create vector of wavenumbers
          [x,y] = meshgrid(0:hx:1,0:hy:1);
          kk    = feval(kvar,x,y);
          kk    = reshape(kk',[nv,1]);
          Ksq   = spdiags(kk,0,nv,nv).*spdiags(kk,0,nv,nv);
          
          %Imaginary shift I*eps
          eps   = feval(epsvar,x,y); size(eps);
          eps   = reshape(eps',[nv,1]);
          Ieps  = 1i*spdiags(eps,0,nv,nv);
          
          %A: 2D Helmholtz matrix
          A = A - Ksq-Ieps;
         
        otherwise
            error('invalid boundary conditions')                  
end

end



