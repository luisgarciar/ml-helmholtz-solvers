%% Helmholtz + CS-MINRES
% In this script we test the performance of the CS-MINRES
% method for the Helmholtz equation vs GMRES

%%  Construction of the Helmholtz Matrix

clc
k    = 50;     %wavenumber
ppw  = 12;      %points per wavelength

nd   = ceil(5*k/pi);    %number of interior discretization points in 1D
nd   = mod(nd,2)*nd + (1-mod(nd,2))*(nd+1); 

flag    = 1;  
dim     = 2;         %dimension    
bc      = 'som';     %type of problem
f       = @rhs;      %right hand side
h       = 1/(nd+1);  %gridsize
[A,ex_sol,b] = helmholtz(f,k,nd,bc,dim,flag); %construction of the 2D matrix
rtol    = 1e-12;
maxit   = 2000;

%% Test of CS-MINRES vs GMRES

tic; 
[sol1,~,~,~,~,~,~,~,~,~,~,resvec1,~] = csminresqlp(A,b,rtol,maxit,[],[],[],[],[],false);
time_csminres = toc;

%semilogy(0:length(resvec1)-1,resvec1/norm(b),'b-.');

tic; 
[sol2,~,~,~,resvec2] = gmres(A,b,[],rtol,maxit);
time_gmres = toc;

semilogy(0:length(resvec1)-1,resvec1/norm(b),'b-.',0:length(resvec2)-1,resvec2/norm(b),'r-.');
hold on


time_gmres
time_csminres


legend('Residuals')
xlabel('Iteration')
ylabel('Relative residuals')


[x,y] = meshgrid(h:h:1-h);

u     = reshape(ex_sol',[nd,nd]);
%gridf = reshape(f',[nd,nd]);

Reu  = real(u);
Imu  = imag(u);

figure
surf(x,y,Reu)


