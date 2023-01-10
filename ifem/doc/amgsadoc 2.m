%% ALGEBRAIC MULTIGRID METHOD: SMOOTHED AGGREGATION
%
% We consider solving an SPD matrix equation |Ax = b|, where |A| could be
% obtained by the finite element discretization on a unstructured grids. A
% coarsening of the graph of |A| is needed and restriction and prolongation
% can be constructued based on the coarsening. We consider smoothe
% aggregation algebraic multigrid based on the strong connectness. A breif
% introduction of sa-AMG can be found at
% <http://link.springer.com/article/10.1007%2FBF02238511 Algebraic
% multigrid by smoothed aggregation for second and fourth order elliptic
% problems>

%% Coarsening
% See <./coarsenAMGadoc.html coarsenAMGadoc> 

%% Prolongation
%
% The prolongation operator is first set to be piecewise constant and then
% using matrix A to smooth out this simple prolongation using weighted
% Jacobi iteration for several steps. The default choice of the weight is
% 0.35 and the smoothingstep is 2. Larger smoothing steps will result a
% denser prolongation operator.

%% Test: Different meshes
%
% We consider linear finite element discretization of the Poisson equation
% with homongenous Dirichlet boundary condition on different meshes. We
% summarize the results in <./amgdoctest1.html AMG Test I> which shows our
% AMG solver is of complexity O(Nlog(logN)).

%% Test: Different boundary conditions
%
% We consider the linear finite element discretization of Poisson equation
% on the unstructured mesh with Dirichlet or Neumann boundary conditions.
% We summarize the results in <./amgdoctest2.html AMG Test II> which shows
% our AMG solver is robust to different boundary conditions.

%% Test: Different time stepsize
%
% We consider the linear finite element discretization of heat equation on
% the unstructured mesh with Neumann boundary conditions. We test the
% implicit time discretization with various time stepsizes dt = 1/h^4,
% 1/h^2, 1/h, 1. We summarize the results in <./amgdoctest3.html AMG Test III> which shows
% our AMG solver is robust to different time steps. 

%% Test: Three dimensional problems
%
% We consider the linear finite element discretization of Poisson equation
% in three dimensions with Dirichlet boundary conditions and compare
% geometric multigrid and algebraic multigrid in <./amgdoctest4.html AMG
% Test IV>.
%
% The algebraic multigrid is robust to the size of the matrix. The
% iteration steps increases slightly. The preformance is better on
% structure grids than unstructure grids. The sparsity of the coarse grid
% is increased
