function [npf,lev] = npcc2npf(npcc,k,ppw)
%% NPCC2F Computes the number of points on the finest grid given the
%         number of points on the coarsest grid and discretization
%         requirements
%
% Use:  [npf,lev] = npcc2npf(npcc,k,ppw)
%
% Input
%       npcc : number of points on coarsest grid
%       k    : wavenumber (for Helmholtz eqn)
%       ppw  : minimum number of points per wavelength   
%
% Output
%       npf : number of points on the finest grid
%       l   : number of levels (grids) from coarsest to finest
%
%  Author:      Luis Garcia Ramos, 
%               Institut fur Mathematik, TU Berlin
%               Version 1.0, Jun 2016
%%  
%The number of grid levels lev and points of the fine grid npf is
%chosen according to the rule (2*pi/k*npc) approx ppw

m   = ppw*k/(2*pi*(npcc+1));
lev = ceil(log2(m));
npf = 2^(lev-1)*(npcc+1)-1;

end

