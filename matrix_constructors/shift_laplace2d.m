function [A, sol,b] = shift_laplace2d(f,kvar,b1,b2,npx,npy,bc,flag)
%% SHIFT_LAPLACE2D: Direct solver for the 2-D  shifted Laplace problem.
%  Solves -div(grad u)-k^2*(b1-i*b2)*u=f with various boundary conditions
%
%  Use: [A, sol,b] =  shift_laplace2d(f,kvar,b1,b2,npx,npy,bc,flag)
%
%  Input: 
%  f:      right-hand side (function handle)
%  k:      wavenumber  of Helmholtz equation (function handle, k depends on (x,y))
%  b1,b2:  real and imaginary shifts 
%  npx:    number of interior discretization points in the x-direction
%  npy:    number of interior discretization points in the y-direction
%  bc:     type of boundary conditions:      
%          'dir' for homogeneous dirichlet bc's
%          'som' for sommerfeld bc's 
%  flag:   if flag==1 solve exactly and return solution
%
%  Output:
%  A:      discrete shifted Laplace operator
%  b:      right hand side vector
%  sol:    solution of the linear system
%
%  Author: Luis Garcia Ramos, 
%          Institut fur Mathematik, TU Berlin
%
%  TO DO:  Add modification for nonuniform grid and nonconstant wavenumber
%
%Version 1.0: Modified for handling nonconstant wavenumbers  Nov 2015
%%%%%

%% Construction of 2D matrix
hx  = 1/(npx+1);            %gridsize
hy  = 1/(npy+1);            %gridsize

switch bc
    case 'dir'
        % Dirichlet 2D matrix (no boundary points)        
        nv = npx*npy;
        v  = ones(nv,1);
        N  = -1/hy^2*v; S=N;
        W  = -1/hx^2*v; E=W;
        C  = (2/hy^2+2/hy^2)*v;
        
        %Create vector of wavenumbers
        [x,y] = meshgrid(hx:hx:1-hx,hy:hy:1-hy);
        kk    = feval(kvar,x,y);
        kk    = reshape(kk',[nv,1]);
        
        %A: 2D Helmholtz matrix
        A = spdiags([S W C E N],[-npx -1 0 1 npx], nv, nv);
        A = A-(b1-1i*b2)*spdiags(kk,0,nv,nv).*spdiags(kk,0,nv,nv);
         
        %b: right hand side constructed with function f 
        b = feval(f,x,y);
        b = reshape(b',[nv,1]);
        
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
             A(n,:)       = zeros(nv,1)';
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
          
          %Adding the effect of the wavenumber k on A
          %Create vector of wavenumbers
          [x,y] = meshgrid(0:hx:1,0:hy:1);
          kk   = feval(kvar,x,y);
          kk   = reshape(kk',[nv,1]);
          bKsq  = (b1-1i*b2)*spdiags(kk,0,nv,nv).*spdiags(kk,0,nv,nv);
          
          %A: 2D Helmholtz matrix
          A = A - bKsq;
          
          %b: right hand side constructed with function f 
          b = feval(f,x,y);
          b = reshape(b',[nv,1]);
         
        otherwise
            error('invalid boundary conditions')                  
end

% solve directly the linear system
sol = zeros(size(b));
if flag == 1
    sol = A\b;  % Solution 
end
