%% Script for testing the function helmholtz1.m

%% Parameters
np   = 200;       %number of interior discretization points in 1D
k    = 0;         %wavenumber
flag = 1;  
dim  = 2;         %dimension    
bc   = 'dir';     %type of problem
f    = @rhs2;     %right hand side
h    = 1/(np+1);  %gridsize

%% Exact solution of Dirichlet problem
% Let u(x,y)=sin(n*pi*x)*sin(m*pi*y)
% Then -divgrad(u)-k^2*u = [(n^2+m^2)*pi^2-k^2]*u=f

% We use this exact solution to check the code
[A1,sol,b1] = helmholtz1(f,k,np,bc,dim,flag);

[x,y] = meshgrid(h:h:1-h);
    u = reshape(sol',[np,np]);
    
 Reu  = real(u);
 Imu  = imag(u);

m=2; n=2; 
exsol = sin(m*pi*x).*sin(n*pi*y); %exact solution of Dirichlet problem

%surf(x,y,Imu-exsol)
h
norm(exsol-Reu)






