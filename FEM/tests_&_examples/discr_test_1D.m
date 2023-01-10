%Test of the discretization of the Helmholtz problem
% -u''- k^2*u-i*eps*u = f in (0,1)
%  u(0) = 0   
%  u'(1)-ku(1)= 0
%
%  References:
% F. Ihlenburg and I. Babuska, Finite element solution of the 
% Helmholtz equation with high wave number Part I: The h-version of the FEM
%
% S Langdon , S Chandler-Wilde, Finite element methods for acoustic scattering 


 

%% Test 1: Fixed wavenumber and grid size, compare exact and discrete sol.
clear all;
k    = 40;
np   = 80;
h    = 1/np; gridx = h*(1:1:np)';
u_ex = exact_sol(k,gridx);
plot(gridx,real(u_ex),'r-');
hold on
A = helmholtzfem(k,np,0);          %Helmholtz matrix
f = ones(np,1); f(np)=0.5; f=h*f;  %Note the scaling in the right hand side!
u = A\f; %u = [0;u];
plot(gridx,real(u),'k-');
relerr = norm(u-u_ex)/norm(u_ex);


%% Test 2: Fixed wavenumber varying gridsize, compare exact and discrete sol.
close all

k        =  10;  %wavenumber
npp      = round(logspace(1,5));
normerr  = zeros(length(npp),1);
relerr   = zeros(length(npp),1);

for j = 1:length(npp)
    np   = npp(j); h=1/np; x=h*(1:1:np)';
    A    = helmholtzfem(k,np,0);         %Helmholtz matrix
    f    = ones(np,1); f(np)=0.5; f=h*f; %Right hand side (scaled)
    M    = mass(np);                     %Mass matrix
    u    = A\f;
    u_ex = exact_sol(k,x);               %Exact (analytic) solution
    err  = u-u_ex;
    normerr(j) = sqrt(real(err'*M*err));
    normu_ex   = sqrt(real(u_ex'*M*u_ex));
    relerr(j)  = normerr(j)/normu_ex;
end
loglog(npp,relerr)
xlabel('Number of Elements')
ylabel('Relative L_2 error (discrete vs exact solution)')
