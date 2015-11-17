%% 1-D Hemholtz and shifted Laplace matrices

%Parameters
np = 800;    %number of interior discretization points
h  = 1/np;  %gridsize
k  = 200;    %wavenumber
b  = 0.5;   %complex shift 

% %% Dirichlet Boundary Conditions
% 
% %1-D Helmholtz with Dirichlet boundary conditions (DBC)
% l    = ones(np,1)*(-1/h^2);  %lower(=upper) diagonal
% d    = ones(np,1)*(2/h^2); 
% A_d  = spdiags([l d l],[-1 0 1],np,np)- k^2*speye(np); 
% 
% %1-D Shifted Laplacian with DBC
% l    = ones(np,1)*(-1/h^2);  %lower(=upper) diagonal
% d    = ones(np,1)*(2/h^2); 
% M_d  = spdiags([l d l],[-1 0 1],np,np)- k^2*(1-1i*0.5)*speye(np); 
% 
% %Preconditioned matrix (with DBC)
% S_d = M_d\A_d;
% eigv= eigs(S, length(S)-2);
% plot(eigv)



%% Sommerfeld boundary conditions
%1-D Helmholtz with Sommerfeld boundary conditions
d     = ones(np,1)*(2-k^2*h^2); 
a     = (1-(k^2*h^2)/2+1i*k*h); 
b     = (1-(k^2*h^2)/2-1i*k*h);
d     =  [a; d; b];   %diagonal
u     =  -ones(np+2,1);      
A_s     = 1/(h^2)*spdiags([u d u],[-1 0 1],np+2,np+2); %Helmholtz matrix

%1-D Shifted Laplacian with Sommerfeld boundary conditions
d     = ones(np,1)*(2-k^2*(1-0.5*1i)*h^2); 
gamma = (1-(k^2*(1-0.5*1i)*h^2)/2-1i*k*h); 
d     =  [gamma; d; gamma];   
u     =  -ones(np+2,1);      
M_s    = 1/(h^2)*spdiags([u d u],[-1 0 1],np+2,np+2); %Shifted Laplace matrix

%Preconditioned matrix
S_s = M_s\A_s;
eigv= eigs(S_s, length(S_s)-2);
plot(eigv,'+'), axis([-1 2 -1 2]), shg
