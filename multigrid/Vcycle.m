function [v] = Vcycle( A, v, f, sigma, nu, mu, w, levs )
%VCYCLE Solves Helmholtz equation using a multigrid V-cycle.
%   %  INPUT:
%        A:     Discrete Helmholtz operator on the finest grid
%        v:     Initial guess
%        f:     right-hand side function
%        sigma: wavenumber
%        nu:    number of presmoothing steps
%        mu:    number of postsmoothing steps
%        levs:  number of multigrid levels

%  OUTPUT       
%        v:     solution computed with the multigrid V-cycle   
    
 
  lmin  = levs;      %number of levels
  nt    = length(f); %length of initial vector
  l     = log2(nt+1);

    if l == lmin %If on coarse level, solve exactly
        [~,v]=helmholtz_1D(f,sigma,1);
    else
        fprintf('Welcome to Level number %d\n',l);
        %Presmoothing and computation of the residual
        v   = wJacobi(A,v,nu,f,w);
        r   = f-A*v;
        
        %Restriction to coarse grid
        fc  = fwrestriction(r); 
        nc  = length(fc);
        vc  = zeros(nc,1);
        
        %Calling Vcycle to solve the error equation
        [Ac, ~] = helmholtz_1D(fc,sigma,0);
         vc     = Vcycle(Ac, vc, fc, sigma, nu, mu, w, levs)
        
        %Prolongation and correction on current grid
         v      = v + lininterpol(vc);
        
        %Postsmoothing
        v      =  wJacobi(A,v,mu,f,w);

    end


end

