function [A] = shift_laplace2(k,b1,b2,npx,npy,bc)
%% SHIFT_LAPLACE2: Constructs matrices for the 2D Helmholtz problem.
%  Constructs the finite difference matrix corresponding
%  to the discretization of the 2D problem u''- k^2(b1-ib2)u=f
%  on the interval (0,1) 
%
%  When homogeneous Dirichlet boundary conditions are used, 
%  the boundary points are eliminated of the linear system.
%  In case of Sommerfeld boundary conditions, the boundary
%  points are included.
%
%  Use: [A] = shift_laplace2(k,b1,b2,npx,npy,bc)
%
%  Input: 
%  k:      wavenumber 
%  b1,b2:  real and imaginary parts of shift
%  np:     number of interior discretization points
%  bc:     type of boundary conditions:      
%          'dir' for homogeneous dirichlet bc's
%          'som' for sommerfeld bc's 
%
%  Output:
%  A:      discretization matrix of shifted Laplace problem
%
%  Author: Luis Garcia Ramos, 
%          Institut fur Mathematik, TU Berlin
%
%  TO DO:  Add modification for nonconstant wavenumber 
%
% Version 0.1 - Nov 2015
%%%%%

%% Construction of 2D matrix

hx  = 1/(npx+1);            %gridsize in x-direction
hy  = 1/(npy+1);            %gridsize in y-direction

switch bc
    case 'dir'
        
        nv = npx*npy;
        v  = ones(nv,1);
        N  = -1/hy^2*v; S=N;
        W  = -1/hx^2*v; E=W;
        C  = (-k^2*(b1+1i*b2)+2/hy^2+2/hy^2)*v;
        
        A = spdiags([S W C E N],[-npx -1 0 1 npx], nv, nv);
       
   case 'som'
        %2D matrix with Sommerfeld bc's (with boundary points)
        nv = (npx+2)*(npy+2);
        
        %Interior points (i*hx,j*hy)
        v = ones(nv,1);
        N = -1/hy^2*v; S=N;
        W = -1/hx^2*v; E=W;
        C = (-k^2*(b1-1i*b2)+2/hy^2+2/hy^2)*v;
        
        A = spdiags([S W C E N],[-(npx+2) -1 0 1 (npx+2)], nv, nv);         
        
        %Noncorner boundary points
        
        %West boundary: (0,j*hy)
         for j = 1:npy
              k = k;
              n = j*(npx+2)+1; %changing from (0,j) to global index
              A(n,:)       = zeros(nv,1)'; %clear row
              A(n,n)       = (-k^2*(b1-1i*b2)-2*1i*k/hx+2/hx^2+2/hy^2);
              A(n,n+1)     = -2/hx^2;
              A(n,n+npx+2) = -1/hy^2;
              A(n,n-npx-2) = -1/hy^2;   
         end
         
             
         %South boundary: (ihx,0)
         for i=1:npx
              k=k;
              n=i+1; %changing from (i,0) to global index
              A(n,:)       = zeros(nv,1)';
              A(n,n)       = (-k^2*(b1-1i*b2)-2*1i*k/hy+2/hx^2+2/hy^2);
              A(n,n+1)     = -1/hx^2;
              A(n,n-1)    =  -1/hx^2;
              A(n,n+npx+2) = -1/hy^2;  
         end
         
          
         %East boundary: (1,j*hy)
         for j=1:npy
             k = k;
             n = (npx+2)*(j+1); %changing from (npx+1,j) to global index
             A(n,:)       = zeros(nv,1)';
             A(n,n)       = (-k^2*(b1-1i*b2)-2*1i*k/hx+2/hx^2+2/hy^2);
             A(n,n-1)     = -1/hx^2;
             A(n,n-npx-2) = -1/hy^2;  
             A(n,n+npx+2) = -1/hy^2;
         end

         %North boundary: (i*hx,1)
         for i=1:npx
             k = k;
             n = (npx+2)*(npy+1)+(i+1);
             A(n,:)       = zeros(nv,1)';
             A(n,n)         = (-k^2*(b1-1i*b2)-2*1i*k/hy+2/hx^2+2/hy^2);
             A(n,n-1)       =  -1/hx^2;
             A(n,n+1)       =  -1/hx^2;  
             A(n,n-npx-2) = -2/hy^2;
         end
         
          %SW Corner: (0,0)
          k = k;
          A(1,:)       = zeros(nv,1)';
          A(1,1)     = (-k^2*(b1-1i*b2)-2*1i*k/hy-2*1i*k/hx+2/hx^2+2/hy^2);
          A(1,2)     = -2/hx^2;
          A(1,npx+3) = -2/hy^2;
          
          %SE Corner: (1,0)
          k = k;
          n = npx+2;
          A(n,:)    = zeros(nv,1)';
          A(n,n)    = -k^2*(b1-1i*b2)-2*1i*k/hy-2*1i*k/hx+2/hx^2+2/hy^2;
          A(n,n-1)  = -2/hx^2;
          A(n,2*n)  = -2/hy^2;
          
          %NW Corner: (0,1)
          k = k;
          n = (npx+2)*(npy+1)+1;
          A(n,:)         = zeros(nv,1)';
          A(n,n)         = -k^2*(b1-1i*b2)-2*1i*k/hy-2*1i*k/hx+2/hx^2+2/hy^2;
          A(n,n+1)       = -2/hx^2;
          A(n,n-(npx+2)) = -2/hy^2;

          
          %NE Corner: (1,1)
          k = k;
          n = (npx+2)*(npy+2);
          A(n,:)         = zeros(nv,1)';
          A(n,n)         = -k^2*(b1-1i*b2)-2*1i*k/hy-2*1i*k/hx+2/hx^2+2/hy^2;
          A(n,n-1)       = -2/hx^2;
          A(n,n-(npx+2)) = -2/hy^2;

        otherwise
            error('invalid boundary conditions')
 
end

end




