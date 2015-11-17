function [v, sol, error ] = twogridcycle( f, sigma, nu, mu )
%TWOGRIDCYCLE Solves a 1-D Helmholtz equation with a two-grid scheme.
%  INPUT:
%        f:     right-hand side function without values on the boundary
%        sigma: wavenumber
%        nu:    number of presmoothing steps
%        mu:    number of postsmoothing steps
%  OUTPUT       
%        v:     exact solution computed with the backslash operator
%        sol:   solution computed with two-grid cycle    
%        error: norm(err-v,inf);
%
%f = -2*pi*exp(x).*cos(pi*x)+(pi^2+sigma-1)*exp(x).*sin(pi*x);


%% Creates the discrete Helmholtz Operator and stores it in matrix A:

[A, v] = helmholtz_1D(f,sigma); % discretization matrix and exact solution
sol    = zeros(length(f)-2,1);     % prealocates memory for the approximate solution

%% Presmoothing (Apply Jacobi nu-times on fine grid):

w      = 0.666;      % damping parameter for wJacobi smoothing
sol    = wJacobi(A,sol,nu,f(2:length(f)-1),w);
rf     = f(2:length(f)-1)-A*sol; % calculate fine grid residual

%% Restrict fine grid residual to coarse grid:

rc        = fwrestriction(rf);
[~, errc] = helmholtz_1D(rc,sigma); %Solve error equation  on coarse grid

%% Prolongate error to fine grid and correct:

errf = lininterpol(errc);
sol  = sol + errf;

%% Postsmoothing (Apply Jacobi mu times on fine grid):

sol   = wJacobi(A,sol,mu,f(2:length(f)-1),w);

%% Apply Boundary conditions:

sol  = [0 sol' 0]'; % boundary conditions for plotting
error = norm(v-sol,inf);


end

