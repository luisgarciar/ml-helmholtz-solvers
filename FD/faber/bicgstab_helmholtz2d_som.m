% Experiments with the CSL-preconditioned 1D Helmholtz equation.
% We compare timing and the iteration numbers for solving
%    A M^{-1}x = b
% with gmres, for different implementations of A*M^{-1}*x:
%   - lu: A*(U\(L\x))
%   - ilu: A*(iU\(iL\x))
%   - multi-grid
%
% Further, we do the same when using polynomial acceleration with a truncated
% Faber series s_n(z) of 1/z on a bw set. Then the system to solve is
% AM^{-1}s_n(AM^{-1}) x = b.
%
% OBSERVATIONS: (Dirichlet boundary conditions)
%   acceleration).
% - With truncated Faber series less iterations are needed.
% - With truncated Faber series less time is needed starting at k = 50.

clear all;
save_flag = 0;  % save_flag=1: save plots and table, =0 do not save.

%% Setup parameters
%Setup list of wavenumbers
wavenum = 120:20:160;

%wavenum = 50; %% run this when testing changes in the code
%number of interior points in coarsest grid in one dim
npc = 1;
bc  = 'som';
ppw = 12;   %number of points per wavelength

warning off;

%Parameters for GMRES
restart = [];
tol     = 1e-8;
maxit   = 300;

numruns =  2; %number of runs for every experiment (to average later)

%memory allocation for time, residuals
num_dofs    = zeros(length(wavenum),1);
time_mg     = zeros(length(wavenum),1);
time_mg_FS  = zeros(length(wavenum),1);
iter_mg     = zeros(length(wavenum),2);
iter_mg_FS  = zeros(length(wavenum),2);

resvec_mg      = zeros(maxit,1);
resvec_mg_FS   = zeros(maxit,1);

for kk = 1:length(wavenum)
    %wavenumber
    k = wavenum(kk); eps = 0.5*k^2;
    
    %number of points in finest grid and number of levels
    [npf,numlev] = fd_npc_to_npf(npc,k,ppw);
    
    %2D
    dim  = 2;
    A = helmholtz2(k,0,npf,npf,bc);
    M = helmholtz2(k,eps,npf,npf,bc);
    op_type = 'gal';
    [mg_mat,mg_split,restrict,interp] = mg_setup(k,eps,op_type,npc,numlev,bc,dim);
    
    b  = zeros(length(M),1); ind = floor(length(M)/2);  b(ind)=1;
    x0 = zeros(size(b));
    
    %multigrid and gmres parameters
    npre = 1; npos = 1; w = 0.6; smo = 'wjac'; numcycles = 1;       
    f_mg  = @(v) A*feval(@Fcycle,mg_mat,mg_split,restrict,interp,x0,v,npre,npos,w,smo,1);
    
    %% Polynomial acceleration with truncated Faber series
    
    deg = 1; % degree of truncated Faber series of 1/z.
   
    %Parameters for the bratwurst shaped set
    lambda    = -1; % so that 0 is not in the bw set.
    phi       = pi/2;  % or phi = 0.1*pi; ??
    eps_thick = 0.005;
    [~, ~, ~, M_bw, N_bw] = bw_map(lambda, phi, eps_thick);
    
    %Setup LU_FS and MG_FS preconditioners
    f_mg_FS  = @(x)truncFS(f_mg,f_mg(x),deg,M_bw,N_bw,'fun');
    
    for t=1:numruns                 
        % with f_mg
        tic
        [~,~,~,iter_mg(kk,:),resvec_mg] = bicgstab(f_mg, b, restart, tol, maxit);
        time_mg(kk,1) = time_mg(kk,1)+ toc;
         
        % with f_mg_FS
        tic
        profile on
        [X_mg_FS,~,~,iter_mg_FS(kk,:),resvec_mg_FS] = bicgstab(f_mg_FS,b,restart,tol,maxit);
        time_mg_FS(kk,1) = time_mg_FS(kk,1)+ toc;
        profile off
    end
    
end % of going through different wavenumbers
time_mg     = time_mg/numruns;
time_mg_FS  = time_mg_FS/numruns;



%% Plot timings and save the plot in tikz (.tex) and .eps formats
figure(1)
plot(wavenum, time_mg, 'k-')
hold on
plot(wavenum, time_mg_FS, 'b-')
hold off
ylabel('Time (s)')
xlabel('Wavenumber')
legend('CSL preconditioner', 'Faber series with CSL preconditioner', 'Location','NorthWest');
FS = 14; % font size
set(gca,'LooseInset',get(gca,'TightInset'))
set(gca,'FontSize',FS)

%Saving the plot in tikz format
%Uses matlab2tikz, be sure to include it in the matlab path
kmin = num2str(min(wavenum));
kmax = num2str(max(wavenum));
dimn = num2str(dim);
degr = num2str(deg);
bdc  = num2str(bc);

plot_time_tex = strcat('exp1_time_vs_k_',dimn,'D_',...
                        bdc,'.tex');

%get path of current .m file
currentpath  = mfilename('fullpath');

%generate file name, for saving in tex_files/figures/new_exp
path  = fullfile(currentpath,'..','..','..','tex_files',...
    'figures','new_exp');

plot_file_tex = fullfile(path,plot_time_tex);

if (save_flag == 1 )
    %Save as tikz figure in .tex file
    matlab2tikz('filename',plot_file_tex,'standalone',true,'extraaxisoptions',...
        ['xlabel style={font=\Large},', ...
        'ylabel style={font=\Large},']);
    
    %Save the plot as .eps
    plot_time_eps = strcat('exp1_time_vs_k_',dimn,'D_',bdc,'.eps');
    plot_file_eps = fullfile(path,plot_time_eps);
    print('-depsc2',plot_file_eps)
end


%% Plot iteration numbers and save as tikz (.tex) and .eps
figure(2)
plot(wavenum, iter_mg(:,2),'k-');
hold on
plot(wavenum, iter_mg_FS(:,2),'b-');
hold off
ylabel('Number of GMRES iterations');
xlabel('Wavenumber');
legend('CSL+MG preconditioner', ...
       'Faber series with CSL+MG preconditioner', 'Location','NorthWest');
FS = 14;
set(gca,'LooseInset',get(gca,'TightInset'))
set(gca,'FontSize',FS)

plot_iter_tex = strcat('exp1_iter_vs_k_',dimn,'D_',...
    bdc,'.tex');


%get path of current .m file
currentpath   = mfilename('fullpath');
%generate file name, for saving in tex_files/figures/new_exp
path          = fullfile(currentpath);
plot_file_tex = fullfile(path,plot_iter_tex);


if (save_flag == 1)
    %Save as tikz figure in .tex file
    matlab2tikz('filename',plot_file_tex,'standalone',true,'extraaxisoptions',...
        ['xlabel style={font=\Large},', ...
        'ylabel style={font=\Large},']);
    
    %Save the plot as .eps
    plot_iter_eps = strcat('exp1_iter_vs_k_',dimn,'D_',bdc,'.eps');
    plot_file_eps = fullfile(path,plot_iter_eps);
    print('-depsc2',plot_file_eps);
end

