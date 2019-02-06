function [fovAhat] = computefovAhat(kmult,opts)
%computefovAhat Computes the field of values of 
%preconditioned Helmholtz matrices

%Input:
%kmult: vector of integers 
        %(wavenumbers, to be multiplied by k)
%opts: Data structure with options for the computation
%opts.prec = 'csl' or 'adef': preconditioner
%opts.dim  = dimension of the problem (1 or 2)
%opts. bc  = Boundary conditions, 'som' or 'dir' 
%opts.disc = Discretization options 
%            'q' for quasioptimal --> n = Ck^2
%            'pf' for pollution free --> n = Ck^1.5
%            'ppw' for points per wavelength --> n = Ck  
%
%opts.fvpts = points in the field of values




%   Detailed explanation goes here
outputArg1 = inputArg1;
outputArg2 = inputArg2;
end

