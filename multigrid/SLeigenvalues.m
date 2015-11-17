clear all
close all

%% Parameter definition
n     = 2^3;    %number of interior gridpoints
w     = 2/3;    %damping parameter for weighted Jacobi smoothing
sigma = 50;      %wavenumber
maxit = 200;    %maximum number of Jacobi iterations
nu    = 2;      %number of presmoothing steps
mu    = 2;      %number of postsmoothing steps
lev   = 3;      %number of multigrid levels

beta1 = 1;      %parameters of Shifted Laplacian Preconditioner
beta2 = 0.5;  


%% Construction of the grid and computation of exact solution
x        = [0:1/n:1]';         %1-D grid
x_res    = x(2:length(x)-1);
f        = zeros(size(x));     %rhs function
f((n/2)+1,1) = 10; 
f        = f(2:length(f)-1,1);       %rhs function restricted to inner grid points

[A, sol] = helmholtz_1D(f,sigma,1); %discretization matrix without boundary conditions and exact solution
%sol      = [0 sol' 0]';            %exact solution with boundary conditions
[S, ~ ] = SLaplace_1D(f,sigma,beta1,beta2,flag);

invS = inv(S);
M    = invS*A;
lambda = eigs(M);


%v = zeros(length(f),1);
 

%% Construction of Restriction and Interpolation Matrices

I=eye(n-1);
E=zeros(n/2-1,n-1);
EH=zeros(n-1,n/2-1);
Id=eye(n/2-1);

for i=1:length(I)
    E(:,i)  =  fwrestriction(I(:,i));
end

for j=1:length(Id)
    EH(:,j) = lininterpol(Id(:,j));
end




