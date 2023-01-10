% Experiments with the CSL-preconditioned 2D Helmholtz equation.
% We compare timing and the iteration numbers for solving
%    AM^{-1}x = b
% The 2D Sommerfeld problem is discretized with first order
% finite differences

clear all;
save_flag = 0;  % save_flag=1: save plots and table, =0 do not save.

%% Setup parameters

%Setup list of wavenumbers
%wavenum = [20:20:120];
%wavenum = [5 10 20 40 60 80 100 120 150];
wavenum =  [10 20 40 80];

bc  = 'som1';
ppw = 0.5;   %pollution free grid
warning off

%Parameters for GMRES
restart = [];
tol     = 1e-10;
maxit   = 100;
numruns =  1; %number of runs for every experiment (to average later)

%memory allocation for time, residuals
num_dofs      = zeros(length(wavenum),1);
time_mg       = zeros(length(wavenum),1);
iter_mg       = zeros(length(wavenum),2);
resvec_mg     = zeros(maxit,1);

profile on

for kk = 1:length(wavenum)
    %wavenumber
    k = wavenum(kk);
    factoreps = 0.5;
    dim       = 2;     %2D
    poweps    = 2;
    eps       = factoreps*k^poweps; 
    npcc      = 4;   %number of points in the coarsest grid in 1D
    ppw       = 12; %ppw < 1 for pollution free grid
    
    %number of points in finest grid and number of levels
    [npf,numlev] = fd_npc_to_npf(npcc,k,ppw);
    
    %Helmholtz matrix in 2D
    A = helmholtz2(k,0,npf,npf,bc);
    
    %Multigrid structure for shifted Laplacian
    op_type = 'gal';
    [mg_mat_SL,mg_split,restr,interp] = mg_setup(k,eps,op_type,npcc,numlev,bc,dim);
    
    %Shifted Laplace matrix
    Aeps = mg_mat_SL{1};
    
    %Right hand side
    b  = zeros(length(A),1); ind = floor(length(b)/2);  b(ind)=1;
    x0 = zeros(size(b));
    
    %multigrid and gmres parameters
    npre  = 1; npos = 1; w = 0.5; smo = 'wjac'; numcycles = 1;
    f_mg  = @(v) A*feval(@Vcycle,mg_mat_SL,mg_split,restr,interp,x0,v,npre,npos,w,smo,1);
      
    for t=1:numruns
        % with f_mg
        tic
        [~,~,~,iter_mg(kk,:),resvec_mg] = gmres(f_mg, b, restart, tol, maxit);
        time_mg(kk,1) = time_mg(kk,1)+ toc;        
    end
    
end % of going through different wavenumbers
time_mg  = time_mg/numruns;
profile off

%% Count matrix-vector operations
% For GMRES on A*x = b, each step requires one multiplication by the matrix A.
% For GMRES on p(A)A*x = b, each step requires one multiplication by the matrix
% p(A)A, i.e., 1 + deg(p) multiplications by A.
% Applying this to the matrix M^{-1} A in our problem, we see that LU, iLU, MG
% have ITER multiplications by M^{-1}A, i.e., M^{-1} is applied ITER times.
% For LU+FS, iLU+FS, MG+FS, we have ITER * (1+d) multiplications by M^{-1}A and
% applications of M^{-1}.

mvop_mg  = iter_mg(:,2);

%% Plot timings 
figure(1)
plot(wavenum, time_mg, 'k-*')
ylabel('Time (s)')
xlabel('Wavenumber')
legend('CSL (MG)','Location','NorthWest')

FS = 16; % font size
set(gca,'LooseInset',get(gca,'TightInset'))
set(gca,'FontSize',FS)
box on

%% Plot iteration numbers 
figure(2)
plot(wavenum, iter_mg(:,2), 'k-*')
ylabel('Number of GMRES iterations')
xlabel('Wavenumber')
legend('CSL (MG)','FontSize',FS,'Location','NorthEast')
%title(['2D Helmholtz with CSL-preconditioner and poly. accel., deg=',...
%num2str(deg)])
FS = 14; % font size
set(gca,'LooseInset',get(gca,'TightInset'))
set(gca,'FontSize',FS);

%% Plot relative GMRES residuals vs number of applications of M^{-1}
% for fixed wavenumber (the last one). For a different wavenumber: grab the
% vector with GMRES residuals and modify the following plot.
%
% figure(3)
% semilogy(0:mvop_lu, resvec_lu/resvec_lu(1), 'k--')
% hold on
% % semilogy(0:mvop_ilu, resvec_ilu/resvec_ilu(1), 'k-.')
% semilogy(0:mvop_mg, resvec_mg/resvec_mg(1), 'k-')
% semilogy(0:(1+deg):mvop_lu_FS, resvec_lu_FS/resvec_lu_FS(1), 'b--')
% % semilogy(0:(1+deg):mvop_ilu_FS, resvec_ilu_FS/resvec_ilu_FS(1), 'b-.')
% semilogy(0:(1+deg):mvop_mg_FS, resvec_mg_FS/resvec_mg_FS(1), 'b-')
% hold off
% ylabel('Relative residual in GMRES')
% xlabel('Number of applications of M^{-1}')
% legend('Exact inversion of preconditioner', ...
%        'MG preconditioner', ...
%     'Faber series with exact preconditioner', ...
%     'Faber series with MG preconditioner', ...
%     'Location','NorthWest')
% %title(['1D Helmholtz with CSL-preconditioner and poly. accel., deg=',...
% %num2str(deg)])
% FS = 14; % font size
% set(gca,'LooseInset',get(gca,'TightInset'))
% set(gca,'FontSize',FS)
