%% Script for testing the functions helmholtz1d.m  

%% Parameters
np   = 200;       %number of interior discretization points in 1D
k    = 50;        %wavenumber
flag = 1;
bc   = 'som';     %type of problem
f    = @dirac1d;  %right hand side
h    = 1/(np+1);  %gridsize

[A,u,b] = helmholtz1d(f,k,np,bc,flag);
    
 Reu  = real(u);
 Imu  = imag(u);
 
% Sommerfeld boundary conditions
 x  = 0:h:1;
 plot(x,Reu)
 plot(eig(full(A)),'+r')
 
 
 
 
 



