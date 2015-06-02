%% Helmholtz + CS-MINRES
% In this script we test the performance of the CS-MINRES
% method for the Helmholtz equation vs GMRES

%%  Construction of the Helmholtz Matrix

k    = 50;     %wavenumber
ppw  = 10;     %points per wavelength

nd   = ceil(5*k/pi);    %number of interior discretization points in 1D
nd   = mod(nd,2)*nd + (1-mod(nd,2))*(nd+1); 

flag    = 0;  
dim     = 2;         %dimension    
bc      = 'som';     %type of problem
f       = @rhs;      %right hand side
h       = 1/(nd+1);  %gridsize
[A,~,b] = helmholtz(f,k,nd,bc,dim,flag); %construction of the 2D matrix
rtol    = 1e-6;
maxit   = 2000;

%% Test of CS-MINRES vs GMRES

[~,~,~,~,~,~,~,~,~,~,~,resvec1,~] = csminresqlp(A,b,rtol,maxit);
[x,~,~,~,resvec2] = gmres(A,b,[],rtol,maxit);

semilogy(0:length(resvec1)-1,resvec1/norm(b),'b-.',0:length(resvec2)-1,resvec2/norm(b),'r-.');
hold on

legend('Residuals')
xlabel('Iteration')
ylabel('Relative residuals')

