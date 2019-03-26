%% Construction of the matrices
close all

dim       = 1;
poweps    = 2;
factoreps = 1;
bc        = 'som';

%Wavenumber
k  =  50;
iter_SL = zeros(length(k),2);
minfov  = zeros(length(k),1);

%Line colors and types for plots
linetyp = {'-','-.',':','--','-'};
%color   = {'r','b','g','k'};
color1  = [0 0 0; 0.5 0 0.5; 0 0 1; 0 0.5 0; 1 0 0];
marker  = {'*','o','.','x'};
opt     = {'m','b','g','k'};
%str2={'k','k.-','k--','b*-'};

plot_csl     = 'yes';
plot_defcsl  = 'yes';
ppw = 12;
eps = factoreps*k^poweps;



npc = ceil(k*ppw/pi);
npf = 2*npc+1;

A    = helmholtzfem(k,npf,0,bc);           %Helmholtz matrix
Aeps = helmholtzfem(k,npf,eps,bc);         %Shifted Laplace matrix
M    = mass(npf,bc);

[L,U] = lu(Aeps);
LH    = U'; UH = L';
N     = length(A);

%Let C = MAepsinv
C    = @(x) M*(U\(L\x));
CH   = @(x) UH\(LH\(M*x));
HerC = @(x) 0.5*(feval(C,x) + feval(CH,x)); %Hermitian part of C


opts.isreal = 0;
opts.p      = 60;
opts.issym  = true;
opts.tol    = 1e-10;
opts.v0     = rand(N,1);
fprintf('beginning computation of max eigvalue of H(C) for k = %d  \n', k);
tic
[vmaxHerC, eigmaxHerC] =  eigs(HerC,N,1,'LM',opts);
time_maxeigv = toc;
fprintf('computation of max eigenvalue of H(C) for k=%d finished \n', k);
fprintf('time maxeigv: %f  \n', time_maxeigv);


%Use the maximum eigenvalue of the Hermitian part of C
%to initiate computation of fovC
fprintf('beginning computation of fov for k=%d \n', k);
tic
[fovC,~,~] = sfov(C,CH,vmaxHerC,N,60);
time_fov = toc;
fprintf('time fov: %f  \n', time_fov);

fovAhat = 1+1i*eps*fovC;
reFOV  = real(fovAhat); imFOV = imag(fovAhat);
cvh    = convhull(reFOV,imFOV);

%plotting the fov
label       = strcat('$k = ', num2str(k),'$');
fovplot  = plot(reFOV,imFOV,'Color','k',...
    'LineWidth',2.5,'linestyle','-',...
    'DisplayName',label);

hold on
plot(0,0,'k.','Markersize',5,'LineWidth',2);
plot(1,0,'b.','Markersize',5,'LineWidth',2);
axis equal
axis([-0.2 1.2 -0.7 0.7]);
xlabel('Re(z)','FontSize',16,'Interpreter','latex');
ylabel('Im(z)','FontSize',16,'Interpreter','latex');
h=gca;
set(gca,'Xtick',[0 0.5 1],'FontSize',16);
set(gca,'Ytick',[-0.5 0 0.5],'FontSize',16);
set(gca,'TickLabelInterpreter', 'tex');

%%


