function [A] = helmholtz2var(kvar,epsvar,npx,npy,bc)
%% HELMHOLTZ2VAR: Constructs matrices for the 2D Helmholtz 
%  and shifted Laplace problems with variable wavenumbers.
%
%  Use: [A] = helmholtz2(kvar,epsvar,npx,npy,bc)
%
%  Constructs the finite difference matrix corresponding
%  to the discretization of 
%
%  div(grad u(x))- (k(x)^2+i*eps(x))u(x) = f in (0,1)x(0,1)
%
% With boundary conditions
%  u  =   0 on boundary  (Dirichlet)
%  or
%  du/dn-iku = 0 (Sommerfeld)
%
% When homogeneous Dirichlet boundary conditions are used, 
% the boundary points are eliminated of the linear system.
% In case of Sommerfeld boundary conditions the boundary points
% are included.
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
% Version 1.0 - Nov 2016
%%%%%
%% Construction of 2D matrix
hx  = 1/(npx+1);            %gridsize
hy  = 1/(npy+1);            %gridsize

switch bc
    case 'dir'
        % Dirichlet 2D matrix (no boundary points)        
        np = npx*npy;
        v  = ones(np,1);
        N  = (-1/hy^2)*v; S=N;
        W  = (-1/hx^2)*v; E=W;
        C  = (2/hx^2+2/hy^2)*v;
        
        %Create vector of wavenumbers
        [x,y] = meshgrid(hx:hx:1-hx,hy:hy:1-hy);
         kk   = feval(kvar,x,y);
         kk   = reshape(kk',[np,1]);
         eps  = feval(epsvar,x,y);
         eps  = reshape(eps',[np,1]);
        
        %A: 2D Helmholtz matrix
        A = spdiags([S W C E N],[-npx -1 0 1 npx], np, np);
        
        for i=1:(npy-1)       %Modify points closest to east and west boundaries
              ii=npx*(i-1)+npx;
              A(ii,ii+1) = 0;
              A(ii+1,ii) = 0;
        end
        
        A = A-(spdiags(kk,0,np,np).*spdiags(kk,0,np,np))...
              -1i*spdiags(eps,0,np,np);
         
    case 'som'        
        %2D matrix with Sommerfeld bc's (with boundary points)
        nx = npx+2;
        ny = npy+2;
        np = nx*ny;
        
        %Interior points (i*hx,j*hy)
        v = ones(np,1);
        N = -1/hy^2*v;  S = N;
        W = -1/hx^2*v;  E = W;
        C = (2/hy^2+2/hy^2)*v;
        
        %A: Negative Laplacian (Helmholtz matrix without the 0-order term)
        A = spdiags([S W C E N],[-(npx+2) -1 0 1 (npx+2)], np, np); 
        
        %Corner points
        %SW Corner: (0,0)
        k = feval(kvar,0,0);
        A(1,:)     =  zeros(np,1)';
        A(1,1)     =  2/hx^2+2/hy^2-2*1i*k/hy-2*1i*k/hx;
        A(1,2)     = -2/hx^2;
        A(1,1+nx)  = -2/hy^2;
      
        %SE Corner: (1,0)
        k = feval(kvar,1,0);
        n=nx;
        A(n,:)    = zeros(np,1)';
        A(n,n)   =  2/hx^2+2/hy^2-2*1i*k/hy-2*1i*k/hx;
        A(n,n-1) = -2/hx^2;
        A(n,2*n) = -2/hy^2;
        
        %NW Corner: (0,1)
        k = feval(kvar,0,1);
        n = (nx)*(ny-1)+1;
        A(n,:)    = zeros(np,1)';
        A(n,n)    = 2/hx^2+2/hy^2-2*1i*k/hy-2*1i*k/hx;
        A(n,n+1)  = -2/hx^2;
        A(n,n-nx) = -2/hy^2;

        %NE Corner: (1,1)
        k = feval(kvar,1,1);
        n = nx*ny;
        A(n,:)     =  zeros(np,1)';
        A(n,n)     =  2/hx^2+2/hy^2-2*1i*k/hy-2*1i*k/hx;
        A(n,n-1)   = -2/hx^2;
        A(n,n-nx)  = -2/hy^2;
                
        %Noncorner boundary points     
        %West boundary: (0,j*hy)        
        j   = 1:npy; 
        ind = j*nx+1;
        k = feval(kvar,0,j*hy);
        Wc  =  2/hx^2+2/hy^2-2*1i*k/hx;
        Ws  = -1/hy^2;
        We  = -2/hx^2;
        Wn  = -1/hy^2;
        WC  = sparse(ind,ind,Wc,np,np); 
        WN  = sparse(ind,ind+nx,Wn,np,np);      
        WS  = sparse(ind,ind-nx,Ws,np,np);
        WE  = sparse(ind,ind+1,We,np,np);
        W  = WC+WN+WS+WE;
        A(ind,:)=W(ind,:);   
             
        %South boundary: (ihx,0)         
        ind = (1:npx)+1;
        k   = feval(kvar,(1:npx)*hx,0);
        Sc  = 2/hx^2+2/hy^2-2*1i*k/hy;
        Sw  = -1/hx^2;
        Se  = -1/hx^2;
        Sn  = -2/hy^2;
        SC  = sparse(ind,ind,Sc,np,np);
        SE  = sparse(ind,ind+1,Se,np,np);
        SW  = sparse(ind,ind-1,Sw,np,np);
        SN  = sparse(ind,ind+nx,Sn,np,np);
        S   = SC+SE+SW+SN;
        A(ind,:)=S(ind,:);
                
        %East boundary
        j  = 1:npy; ind = nx*(j+1);
        k  = feval(kvar,1,j*hy);
        Ec = -2*1i*k/hx+2/hx^2+2/hy^2;
        Ew = -2/hx^2;
        Es = -1/hy^2;
        En = -1/hy^2;
        EC = sparse(ind,ind,Ec,np,np);
        EW = sparse(ind,ind-1,Ew,np,np);
        ES = sparse(ind,ind-nx,Es,np,np);
        EN = sparse(ind,ind+nx,En,np,np);
        E  = EC+EW+EN+ES;
        A(ind,:)=E(ind,:);
                
        %North boundary
        i  = 1:npx; ind = (nx)*(npy+1)+(i+1);
        k  = feval(kvar,i*hx,1);       
        Nc = 2/hx^2+2/hy^2-2*1i*k/hy;
        Ne = -1/hx^2;
        Nw = -1/hx^2; 
        Ns = -2/hy^2;
        NC = sparse(ind,ind,Nc,np,np);
        NW = sparse(ind,ind-1,Nw,np,np);
        NS = sparse(ind,ind-nx,Ns,np,np);
        NE = sparse(ind,ind+1,Ne,np,np);
        N  = NC+NW+NS+NE;
        A(ind,:) = N(ind,:);
          
        %Adding the effect of the wavenumber k and the imaginary shift on A
        %Create vector of wavenumbers
        [x,y] = meshgrid(0:hx:1,0:hy:1);
        kk    = feval(kvar,x,y);
        kk    = reshape(kk',[np,1]);
        Ksq   = spdiags(kk,0,np,np).*spdiags(kk,0,np,np);
          
        %Imaginary shift I*eps
        eps   = feval(epsvar,x,y); size(eps);
        eps   = reshape(eps',[np,1]);
        Ieps  = 1i*spdiags(eps,0,np,np);
          
        %A: 2D Helmholtz matrix
        A = A - Ksq-Ieps;
         
        otherwise
            error('invalid boundary conditions')                  
end

end



