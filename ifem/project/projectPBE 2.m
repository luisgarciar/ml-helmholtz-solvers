%% Project: Nonlinear Poisson Boltzmann Equation
%
% The purpose of this project is to implement Newton's method and FAS for
% solving the nonlinear elliptic equation. The example is the nonlinear
% Poisson-Boltzmann equation for the potential u corresponding to a given
% charge density $\rho (x)$ reads
%
% $$-\Delta u + k^2 \sinh (u) = \rho (x)$$
%
% for $x\in \Omega$, and $u |_{\partial \Omega} = g$.
% 
% For $k = 1$ and $\rho = 0$, an exact solution in 1-d is given by 
%
% $$\bar u(s) = \ln \left ( \frac{1+\cos (s)}{1-\cos (s)}\right). $$ 
%
% We consider a 2-d problem on the unit square $\Omega = (0,1)^2$. Let
% $a=(1.0,2.0)/\sqrt{5}$. We choose $k =1$, $\rho$, and $g$ such that
% the exact solution is $u(x) = \bar u(0.1+(x,a)).$

%% Step 1: Linearied Poisson Boltzmann Equation
%
% * Given a current approximation of u, derive the linearized Poisson
% Boltzmann equation (LPBE) at u.
%
% * Assemble the matrix equation for the LPBE. Besides the matrix of
% Laplacian operator, you need to compute the mass matrix corresponding to
% the L2 inner product. You can use three vertices quadrature rule i.e.
%
% $$\int _{\tau} f(x) dx = \frac{1}{3}\sum _{i=1}^3f(x_i)|\tau|.$$ 
%
% Then the mass matrix becomes diagonal. This is known as mass lumping.
%
% * Use direct solver to solve the matrix equation.
%
% * Use multigrid solver (e.g. amg) to solve the matrix equation. You can
% use your own multigrid methods or call amg in ifem.

%% Step 2: Newton's method on uniform grids
% 
% * Implement the Newton's method. Control the relative error of the
% residual in the stopping criteria.
%
% * Change the tolerance or max iteration steps in multigrid solver and
% collect a table of total iteration steps and cpu time for different
% choices of inner iteration.
%
% * Uniform refine the grid and list the iteration steps for different h

%% Step 3: Nonlinear Multigrid: FAS
%
% * Implement the nonlinear Gauss-Seidel smoother.
%
% * Test two level version of FAS.
%
% * Change two level FAS to V-cycle FAS by recrusion.
%
% * Compare the convergence of FAS with Newton's method.