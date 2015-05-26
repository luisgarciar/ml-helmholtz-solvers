%% Script for testing the function helmholtz.m

%% Parameters
nd   = 101;    %number of interior discretization points in 1D
k    = 50;     %wavenumber
flag = 1;  
dim  = 2;      %dimension    
bc   = 'dir';  %type of problem
f    = @rhs;   %right hand side
h    = 1/(nd+1); %gridsize

%% Solution of the Helmholtz equation and postprocessing the solution
[A, sol,b] = helmholtz(f,k,nd,bc,dim,flag);

[x,y] = meshgrid(h:h:1-h)

u     = reshape(sol',[nd,nd]);
%gridf = reshape(f',[nd,nd]);

Reu  = real(u);
Imu  = imag(u);

figure
surf(x,y,Reu)



















