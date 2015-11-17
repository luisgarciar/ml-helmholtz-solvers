function [ v ] = FMG( A, f, sigma, nu, mu, w, levs)
%FMG Solves the 1-D Helmholtz equation using the FMG-Vcycle method.
%   
%  INPUT:
%        A:     Discrete Helmholtz operator on the finest grid
%        v:     Initial guess
%        f:     right-hand side function
%        sigma: wavenumber
%        nu:    number of presmoothing steps
%        mu:    number of postsmoothing steps
%        levs:  number of multigrid levels

%  OUTPUT       
%        v:     solution computed with the FMG V-cycle 


  lmin  = levs;      %number of levels
  nt    = length(f); %length of initial vector
  l     = log2(nt+1);

    if l == lmin %If on coarse level, solve exactly
        v=zeros(nt,1);
    else
        fc=fwrestriction(r);




end

