clear all
close all

%% Parameter definition
n     = 2^10;    %number of interior gridpoints
w     = 2/3;    %damping parameter for weighted Jacobi smoothing
k     = 0;      %wavenumber
maxit = 200;    %maximum number of Jacobi iterations
nu    = 2;      %number of presmoothing steps
mu    = 2;      %number of postsmoothing steps
lev   = 3;      %number of multigrid levels

%% Construction of the grid and computation of exact solution
x        = [0:1/n:1]'; %1-D grid
x_res    = x(2:length(x)-1);
f        = zeros(size(x));%rhs function
f(1,(n/2)+1) = 10; 
f        = f(2:length(f)-1);       %rhs function restricted to inner grid points

[A, sol] = helmholtz_1D(f,k,1); %discretization matrix without boundary conditions and exact solution
sol      = [0 sol' 0]';            %exact solution with boundary conditions

v = zeros(length(f),1);

%% Test of weighted Jacobi function
%v = wJacobi(A,v,200,f(2:length(f)-1),w);
%v = [0 v' 0]'; %solution computed with weighted Jacobi
%xc = [0:2/n:1];            %coarse grid  
%xc = xc(2:1:length(xc)-1);

 
%%  Two-grid scheme
% 
% %Presmoothing: Relax nu times on fine grid
% v=wJacobi(A,v,nu,f,w);
% rf=f-A*v; %fine grid residual
% 
% %Restriction
% rc=fwrestriction(rf);
% [Acoarse, errc] = helmholtz_1D(rc,sigma); %Solve error equation  on coarse grid
% 
% %Prolongate error to fine grid and correct
% errf=lininterpol(errc);
% v=v+errf;
% 
% %Postsmoothing: Relax mu times on fine grid
% v=wJacobi(A,v,mu,f,w);
% 
% v = [0 v' 0]'; %boundary conditions for plotting
% error=norm(v-sol,inf);



 
 %% V-cycle
 [v] = Vcycle(A, v, f, 0, nu, mu, w, lev);
  v  = [0 v' 0]'; % adding boundary conditions for plotting

 
%% Plot of exact solution vs MG solution
plot(x,sol);
hold on
plot(x,v,'r*');

%% Plot of error
%plot(x,norm(v-sol))


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




