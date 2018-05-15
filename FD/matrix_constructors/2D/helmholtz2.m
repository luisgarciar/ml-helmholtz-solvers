function [A] = helmholtz2(k,eps,npx,npy,bc)
%% HELMHOLTZ2: Constructs matrices for the 2D Helmholtz and shifted Laplace problems.
%  Constructs the finite difference matrix corresponding
%  to the discretization of 
%  div(grad u)- (k^2+i*eps)u = f in (0,1)x(0,1)
%
% With boundary conditions
%  u  =   0 on boundary  (Dirichlet)
%  or
%  du/dn-iku = 0 (Sommerfeld)
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

%npx
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
        nx = npx+2;
        ny = npy+2;
        np = nx*ny;
        
        %Interior points
        W = -ones(np,1)/hx^2;  E=W; %Dxx
        N = -ones(np,1)/hy^2;  S=N; %Dyy
        C = ((2/hx^2) + (2/hy^2)-(k^2+1i*eps))*ones(np,1);
        A = spdiags([S W C E N],[-nx -1 0 1 nx],np, np);
          
        %Corner points
        %SW Corner: (0,0)
        A(1,:)     = zeros(np,1)';
        A(1,1)     = 2/hx^2+2/hy^2-k^2-1i*eps-2*1i*k/hy-2*1i*k/hx;
        A(1,2)     = -2/hx^2;
        A(1,1+nx)  = -2/hy^2;
      
        %SE Corner: (1,0)
        n=nx;
        A(n,:)    = zeros(np,1)';
        A(n,n)   =  2/hx^2+2/hy^2-k^2-1i*eps-2*1i*k/hy-2*1i*k/hx;
        A(n,n-1) = -2/hx^2;
        A(n,2*n) = -2/hy^2;

        %NW Corner: (0,1)
        n = (nx)*(ny-1)+1;
        A(n,:)    = zeros(np,1)';
        A(n,n)    = 2/hx^2+2/hy^2-k^2-1i*eps-2*1i*k/hy-2*1i*k/hx;
        A(n,n+1)  = -2/hx^2;
        A(n,n-nx) = -2/hy^2;

        %NE Corner: (1,1)
        n = nx*ny;
        A(n,:)     =  zeros(np,1)';
        A(n,n)     =  2/hx^2+2/hy^2-k^2-1i*eps-2*1i*k/hy-2*1i*k/hx;
        A(n,n-1)   = -2/hx^2;
        A(n,n-nx)  = -2/hy^2;
     
        %Noncorner boundary points
        
        % Vectorized version
        %West boundary
        j   = 1:npy; 
        ind = j*nx+1;
        Wc  = 2/hx^2+2/hy^2-k^2-1i*eps-2*1i*k/hx;
        Ws  = -1/hy^2;
        We  = -2/hx^2;
        Wn  = -1/hy^2;
        WC  = sparse(ind,ind,Wc,np,np); 
        WN  = sparse(ind,ind+nx,Wn,np,np);      
        WS  = sparse(ind,ind-nx,Ws,np,np);
        WE  = sparse(ind,ind+1,We,np,np);
        W  = WC+WN+WS+WE;
        A(ind,:)=W(ind,:);      
      
        %South boundary
        ind = (1:npx)+1;
        Sc  =  2/hx^2+2/hy^2-k^2-1i*eps-2*1i*k/hy;
        Sw  = -1/hx^2;
        Se  = -1/hx^2;
        Sn  = -2/hy^2;
        SC  = sparse(ind,ind,Sc,np,np);
        SE  = sparse(ind,ind+1,Se,np,np);
        SW  = sparse(ind,ind-1,Sw,np,np);
        SN  = sparse(ind,ind+nx,Sn,np,np);
        S   = SC+SE+SW+SN;
        A(ind,:) = S(ind,:);        
%          
        %East boundary
        j=1:npy; ind = nx*(j+1);
        Ec = (-k^2-1i*eps-2*1i*k/hx+2/hx^2+2/hy^2);
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
        i = 1:npx; ind = (nx)*(npy+1)+(i+1);
        Nc = 2/hx^2+2/hy^2-k^2-1i*eps-2*1i*k/hy;
        Ne = -1/hx^2;
        Nw = -1/hx^2; 
        Ns = -2/hy^2;
        NC = sparse(ind,ind,Nc,np,np);
        NW = sparse(ind,ind-1,Nw,np,np);
        NS = sparse(ind,ind-nx,Ns,np,np);
        NE = sparse(ind,ind+1,Ne,np,np);
        N  = NC+NW+NS+NE;
        A(ind,:)=N(ind,:);
        
    case 'som1'
       %Construction of the 2D matrix (version 2)
        %2D matrix with Sommerfeld bc's (with boundary points)
        nx = npx+2;
        ny = npy+2;
        np = nx*ny;
        
        %Interior points
        W = -ones(np,1)/hx^2;  E = W; %Dxx
        N = -ones(np,1)/hy^2;  S = N; %Dyy
        C = ((2/hx^2) + (2/hy^2)-(k^2+1i*eps))*ones(np,1);
        A = spdiags([S W C E N],[-nx -1 0 1 nx],np, np);
          
        %Corner points
        %SW Corner: (0,0)
        A(1,:)     = zeros(np,1)';
        A(1,1)     = 1/hx^2+1/hy^2-k^2-1i*eps-1i*k/hy-1i*k/hx;
        A(1,2)     = -1/hx^2;
        A(1,1+nx)  = -1/hy^2;
      
        %SE Corner: (1,0)
        n = nx;
        A(n,:)   = zeros(np,1)';
        A(n,n)   =  1/hx^2+1/hy^2-k^2-1i*eps-1i*k/hy-1i*k/hx;
        A(n,n-1) = -1/hx^2;
        A(n,2*n) = -1/hy^2;

        %NW Corner: (0,1)
        n = (nx)*(ny-1)+1;
        A(n,:)    =   zeros(np,1)';
        A(n,n)    =   1/hx^2+1/hy^2-k^2-1i*eps-1i*k/hy-1i*k/hx;
        A(n,n+1)  =  -1/hx^2;
        A(n,n-nx) =  -1/hy^2;

        %NE Corner: (1,1)
        n = nx*ny;
        A(n,:)     =  zeros(np,1)';
        A(n,n)     =  1/hx^2+1/hy^2-k^2-1i*eps-1i*k/hy-1i*k/hx;
        A(n,n-1)   = -1/hx^2;
        A(n,n-nx)  = -1/hy^2;
     
        %Noncorner boundary points    
        %Vectorized version
        %West boundary
        j   =  1:npy; 
        ind =  j*nx+1;
        Wc  =  -k^2 - 1i*eps + 2/hy^2 + 1/hx^2 - 1i*k/hx;
        Ws  = -1/hy^2;
        Wn  = -1/hy^2;
        We  = -1/hx^2;
        WC  = sparse(ind,ind,Wc,np,np); 
        WN  = sparse(ind,ind+nx,Wn,np,np);      
        WS  = sparse(ind,ind-nx,Ws,np,np);
        WE  = sparse(ind,ind+1,We,np,np);
        W   = WC+WN+WS+WE;
        A(ind,:)= W(ind,:);      
      
        %South boundary
        ind = (1:npx) + 1;
        Sc  = -k^2-1i*eps + 2/hx^2 + 1/hy^2 -1i*k/hy;
        Sw  = -1/hx^2;
        Se  = -1/hx^2;
        Sn  = -1/hy^2;
        SC  = sparse(ind,ind,Sc,np,np);
        SE  = sparse(ind,ind+1,Se,np,np);
        SW  = sparse(ind,ind-1,Sw,np,np);
        SN  = sparse(ind,ind+nx,Sn,np,np);
        S   = SC+SE+SW+SN;
        A(ind,:) = S(ind,:);        
          
        %East boundary
        j  = 1:npy; ind = nx*(j+1);
        Ec = (-k^2-1i*eps +2/hy^2 +1/hx^2 -1i*k/hx);
        Ew = -1/hx^2;
        Es = -1/hy^2;
        En = -1/hy^2;
        EC = sparse(ind,ind,Ec,np,np);
        EW = sparse(ind,ind-1,Ew,np,np);
        ES = sparse(ind,ind-nx,Es,np,np);
        EN = sparse(ind,ind+nx,En,np,np);
        E  = EC+EW+EN+ES;
        A(ind,:)= E(ind,:);
                      
        %North boundary
        i  = 1:npx; ind = (nx)*(npy+1)+(i+1);
        Nc = - k^2 -1i*eps+ 2/hx^2 + 1/hy^2 -1i*k/hy;
        Ne = -1/hx^2;
        Nw = -1/hx^2; 
        Ns = -1/hy^2;
        NC = sparse(ind,ind,Nc,np,np);
        NE = sparse(ind,ind+1,Ne,np,np);
        NW = sparse(ind,ind-1,Nw,np,np);
        NS = sparse(ind,ind-nx,Ns,np,np);
        N  = NC+NW+NS+NE;
        A(ind,:)=N(ind,:); 
                 
    otherwise
        error('invalid boundary conditions')
end

end



