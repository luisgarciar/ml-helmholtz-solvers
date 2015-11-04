%% Script for testing the function helmholtz1.m

%% Parameters
np   = 150;       %number of interior discretization points in 1D
k    = 0;         %wavenumber
flag = 1;  
dim  = 2;         %dimension    
bc   = 'dir';     %type of problem
f    = @rhs3;     %right hand side
h    = 1/(np+1);  %gridsize

%% Exact solution of Dirichlet problem
% Let u(x,y) given by solrhs3
% Then -divgrad(u)-k^2*u = rhs3

% We use this exact solution to check the code
[A1,sol,b1] = helmholtz1(f,k,np,bc,dim,flag);

[x,y] = meshgrid(h:h:1-h);
    u = reshape(sol',[np,np]);
    
 Reu  = real(u);
 Imu  = imag(u);

m=2; n=2; 
exsol = solrhs3(x,y); %exact solution of Dirichlet problem

surf(x,y,u)
h
norm(exsol+Reu)



