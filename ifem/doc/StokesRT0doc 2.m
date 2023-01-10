%% Stokes Elements in 2-D: RT0-P0 
%
% We explain stable element RT0-P0 (the lowest order Raviart-Thomas element
% for velocity and piecewise constant for pressure) of Stokes equations in
% 2-D.
%
% [u,p,w,edge,eqn,info] = StokesRT0(node,elem,bdFlag,pde,option)
% 
% See also PoissonRT0doc
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

%% Data Structure
%
% The data structure for RT0 element can be found at
% <dofRT0doc.html RT0 lowest order edge element in 2D>.
%
% For pressure, the basis is the chacteristic function of each triangle. So
% in the computation of divergence operator, elemSign should be used. 

[node,elem] = squaremesh([0 1 0 1], 0.5);
[elem,bdFlag] = sortelem(elem,bdFlag);
[elem2edge,edge] = dofedge(elem);

%% Discretization
%
%  We discretize the following varational form of Stokes equation
%  ------------------------------------------------------------------
% (rot u,rot v) +(div u, div v)  - (p, div v)  = (f,v)
%                                 (- div u,q)   = 0   
%                                           u   = g   on boundary
%  ------------------------------------------------------------------
