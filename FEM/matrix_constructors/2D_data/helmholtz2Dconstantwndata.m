function pde = helmholtz2Dconstantwndata(k,factoreps,poweps)
%% HELMHOLTZ2CONSTANTWNDATA 
%  Data for Helmholtz/shifted Laplace problem with 
%  constant wavenumber
%
% Returns a struct with data related to the
% the Helmholtz/shifted Laplace problem
%   
%  -div(grad u)-(k^2 + i*eps)u = 0 in Omega= (0,1)x(0,1)
%         grad(u) dot n - i*ku = g in bd(Omega) 
%
%  Usage:  
%  Input: wavenumber k, factoreps, poweps 
%       (for shifted Laplacian eps= factoreps*k^poweps)
%
%   Output: 
%       pde: struct containing the following data:
%            wavenumber k, factoreps, poweps 
%               
% Created by Luis Garcia Ramos, Aug 2017, based
% on iFEM test files
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.


kk  = k; 
pde = struct('k',kk,'factoreps',factoreps,'poweps',poweps);

end