%%Script for plotting the field of values of SLP+Helmholtz matrices
%%and computing quantities related to GMRES convergence

%% Parameters
np   = 500;       %number of interior discretization points in 1D
k    = 150;       %wavenumber
%flag = 0;  

dim  = 2;         %dimension    
bc   = 'dir';     %type of boundary conditions ('dir': dirichlet, 'som': sommerfeld)
f    = @rhs2;     %right hand side
h    = 1/(np+1);  %gridsize
b1 = 1; b2 = 0.5;    %parameters of the shifted Laplacian preconditioner

if (k*h > pi/5)
    error('grid resolution too low for wavenumber');
    %np = ceil(k*5/pi);
end

%% Computation of matrices
e = ones(np,1);
Lap1d = spdiags([e -2*e e], -1:1, np, np)/h^2;
Helm1d = -Lap1d-k^2*eye(np);            %Helmholtz matrix
ShifL1d = -Lap1d-k^2*(b1-1i*b2)*eye(np); %Shifted Laplacian

M = ShifL1d\Helm1d;    %SLP preconditioned Helmholtz matrix
L = eig(M); %eigenvalues of preconditioned matrix

%% Computation of the eigenvalues using exact formulas
% l      = linspace(1,np,np); 
% eigv   = (2+2*cos(l*pi*h));                         %eigenvalues of 1dlaplacian 
% slpeig = (eigv-k^2*h^2)./(eigv-k^2*h^2*(b1-1i*b2)); %eigenvalues of 1d-slp+helmholtz          

%% Plots of eigenvalues from matrix computation and exact formulas
% %figure(1); clf;
% %plot(slpeig,'g.');
% axis([-.1 2 -1.5 1.5])
% hold on
% %plot(L,'g.');

%pause
%% Computation of the field of values of the preconditioned matrix

figure(1)
[inner,outer]=fov(M);
                             
%Vector showing the FOV bound for GMRES (preconditioned system)
c = log(1-inner*outer);
resvec_bd = (1:1:length(M))'; resvec_bd = (c/2)*resvec_bd;
resvec_bd = exp(resvec_bd);


%% GMRES runs to check if the convergence bound describes the actual behavior
rhs = 4*rand(length(M),1); rhs = rhs/norm(rhs);

% Run GMRES on M, b
[~,~,~,~,resvec1]  = gmres(M,rhs,[],1e-12,np);
[~,~,~,~,resvec2]  = gmres(Helm1d,rhs,[],1e-12,np);
echo on
size(resvec_bd)
resvec_bd = resvec_bd(1:length(resvec1),1);
echo off

resvec_bd = [1 ; resvec_bd];
%size(resvec_bd)

%% Plot results
figure(2)
semilogy(0:length(resvec1)-1,resvec1,'b-.',...
0:length(resvec2)-1,resvec2,'r-.',0:length(resvec_bd)-1,resvec_bd,'gs-' );

legend('Residuals')
xlabel('Iteration')
ylabel('Relative residuals')
title('1D-Helmholtz equation')

%Conclusion of experiment so far: The fov bound does not describe gmres
%convergence...

% Continuation of 1-D Experiment: The deflated case











%% Experiments with 2-D Helmholtz
Lap2d = kron(Lap1d, speye(length(Lap1d))) + kron(Lap1d, speye(length(Lap1d))) ;



