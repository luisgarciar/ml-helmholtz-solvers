function [npf,lev] = npc_to_npf(npcc,k,ppw)
%% NPC_TO_NPF Computes the number of points on the finest grid given the
%         number of points on the coarsest grid and discretization
%         requirements
%
% Use:  [npf,lev] = npc_to_npf(npcc,k,ppw)
%
% Input
%       npcc : number of points on coarsest grid (including endpoint 0)
%       k    : wavenumber (for Helmholtz eqn)
%       ppw  : minimum number of points per wavelength   
%
% Output
%       npf : number of points on the finest grid (including endpoint 0)
%       l   : number of levels (grids) from coarsest to finest
%
%  Author:      Luis Garcia Ramos, 
%               Institut fur Mathematik, TU Berlin
%               Version 1.0, Jun 2016

%%  

%The number of grid levels lev and points of the fine grid npf is
%chosen according to the rule (2*pi/k*npf) approx ppw

m   = ppw*k/(2*pi*npcc);
lev = ceil(log2(m));
lev = max(lev,1);   %at least 1 level
npf = npcc*(2^lev);

end

