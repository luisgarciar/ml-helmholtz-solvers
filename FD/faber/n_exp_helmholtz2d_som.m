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
%wavenum = [10:10:90, 100:100:1500];
wavenum = [20:20:140];
%wavenum = [25];

%wavenum = 50; %% run this when testing changes in the code
%number of interior points in coarsest grid in one dim
npc = 1;
bc = 'som';
ppw = 12;   %number of points per wavelength

warning off

%Parameters for GMRES
restart = [];
tol     = 1e-8;
maxit   = 300;

numruns =  1; %number of runs for every experiment (to average later)

%memory allocation for time, residuals
num_dofs    = zeros(length(wavenum),1);
%time_lu     = zeros(length(wavenum),1);
%time_lu_FS  = zeros(length(wavenum),1);
time_mg     = zeros(length(wavenum),1);
time_mg_FS  = zeros(length(wavenum),1);

%iter_lu     = zeros(length(wavenum),2);
%iter_lu_FS  = zeros(length(wavenum),2);
iter_mg     = zeros(length(wavenum),2);
iter_mg_FS  = zeros(length(wavenum),2);

%resvec_lu      = zeros(maxit,1);
resvec_mg      = zeros(maxit,1);
%resvec_lu_FS   = zeros(maxit,1);
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
    
    %tic
    %[L, U] = lu(M); % LU = M
    %time_setup_lu = toc;
    
    %f_lu  = @(x) A*(U\(L\x));
    f_mg  = @(v) A*feval(@Fcycle,mg_mat,mg_split,restrict,interp,x0,v,npre,npos,w,smo,1);
    
    %% Polynomial acceleration with truncated Faber series
    
    deg = 1; % degree of truncated Faber series of 1/z.
    %Parameters for the bratwurst shaped set
    lambda    = -1; % so that 0 is not in the bw set.
    phi       = pi/2;  % or phi = 0.1*pi; ??
    eps_thick = 0.005;
    [~, ~, ~, M_bw, N_bw] = bw_map(lambda, phi, eps_thick);
    
    %Setup LU_FS and MG_FS preconditioners
    %f_lu_FS  = @(x)truncFS(f_lu,f_lu(x),deg,M_bw,N_bw,'fun');
    f_mg_FS  = @(x)truncFS(f_mg,f_mg(x),deg,M_bw,N_bw,'fun');
    
    for t=1:numruns
%         % with f_lu
%         tic
%         [~, ~, ~, iter_lu(kk,:), resvec_lu] = gmres(f_lu, b, restart, tol, maxit);
%         time_lu(kk,1) = time_lu(kk,1)+toc;
%                 
        % with f_mg
        tic
        [~,~,~,iter_mg(kk,:),resvec_mg] = gmres(f_mg, b, restart, tol, maxit);
        time_mg(kk,1) = time_mg(kk,1)+ toc;
        
        % with f_lu_FS
%         tic
%         [~,~,~,iter_lu_FS(kk,:),resvec_lu_FS] = gmres(f_lu_FS,b,restart,tol,maxit);
%         time_lu_FS(kk,1) = time_lu_FS(kk,1)+ toc;
%                  
        % with f_mg_FS
        tic
        profile on
        [X_mg_FS,~,~,iter_mg_FS(kk,:),resvec_mg_FS] = gmres(f_mg_FS,b,restart,tol,maxit);
        time_mg_FS(kk,1) = time_mg_FS(kk,1)+ toc;
        profile off
    end
    
end % of going through different wavenumbers
%time_lu     = time_lu/numruns  + time_setup_lu;
time_mg     = time_mg/numruns;
%time_lu_FS  = time_lu_FS/numruns + time_setup_lu;
time_mg_FS  = time_mg_FS/numruns;

%% Plot timings and save the plot in tikz (.tex) and .eps formats
figure(1)
plot(wavenum, time_lu, 'k--')
hold on
plot(wavenum, time_mg, 'k-')
plot(wavenum, time_lu_FS, 'b--')
plot(wavenum, time_mg_FS, 'b-')
hold off
ylabel('Time (s)')
xlabel('Wavenumber')
legend('Exact inversion of preconditioner', ...
    'MG preconditioner', ...
    'Faber series with exact preconditioner', ...
    'Faber series with MG preconditioner', 'Location','NorthWest')
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
plot(wavenum, iter_lu(:,2),'k+-');
hold on
plot(wavenum, iter_mg(:,2),'k-');
plot(wavenum, iter_lu_FS(:,2),'b+-');
plot(wavenum, iter_mg_FS(:,2),'b-');
hold off
ylabel('Number of GMRES iterations');
xlabel('Wavenumber');
legend('CSL + LU preconditioner','CSL+MG preconditioner', ...
       'Faber series with CSL+LU preconditioner', ...
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

%% Export .tex table with the timings for several preconditioning methods
% table_data = [wavenum.',iter_lu(:,2),time_lu,iter_lu_FS(:,2),time_lu_FS,...
%     iter_ilu(:,2), time_ilu, iter_ilu_FS(:,2), time_ilu_FS, ...
%     iter_mg(:,2),time_mg, iter_mg_FS(:,2),time_mg_FS];
% 
% numrows = size(table_data,1);
% numcols = size(table_data,2);
% tableCaption = strcat('Comparison of number of GMRES iterations and timings for 2D Sommerfeld problem');
% tableLabel   = strcat(dimn,'D_',bdc);
% 
% %dataFormat: number of decimals for the variables of each column
% %.0f no decimals, .2f 2 decimals, .3e scientific notation
% TF = '%.3f';    % Time Format.
% dataFormat = {'%.0f','%.0f',TF,'%.0f',TF,'%.0f',TF,'%.0f',TF,'%.0f',TF,'%.0f',TF};
% 
% header = ['\begin{tabular}','{',repmat(('c'),1,numcols),'}'];
% table  = {'\begin{table}[h]';'\resizebox{\textwidth}{!}{';...
%     '\centering';header;'\toprule'};
% 
% row1   = {'$k$ & \multicolumn{2}{c}{CSL+LU} & \multicolumn{2}{c}{FS+CSL+LU} & ',...
%     '\multicolumn{2}{c}{CSL+iLU} & \multicolumn{2}{c}{FS+CSL+iLU} & ',...
%     '\multicolumn{2}{c}{CSL+MG} & \multicolumn{2}{c}{FS+CSL+MG}\\'};
% row2  = {'& Iter. & Time(s) &  Iter. & Time(s) &  Iter. & Time(s) &  Iter. & Time(s) &  Iter. & Time(s) &  Iter. & Time(s) \\ \midrule'};
% table = [table;row1';row2'];
% 
% for i=1:numrows
%     for j =1:numcols
%         dataValue = num2str(table_data(i,j),dataFormat{j});
%         if j==1
%             rowStr = dataValue;
%         else
%             rowStr = [rowStr,' & ',dataValue];
%         end
%     end
%     table(end+1) = {[rowStr,' \\']};
% end
% 
% footer = {'\end{tabular}}';['\caption{',tableCaption,'}']; ...
%     ['\label{table:',tableLabel,'}'];'\end{table}'};
% table = [table;'\bottomrule';footer];
% 
% % Save the table to a file in folder ../../../tex_files/tables
% namefile = strcat('table_',dimn,'D_',bdc,'.tex');
% currentpath = mfilename('fullpath');
% f = fullfile(currentpath,'..','..','..','tex_files','tables',namefile);
% 
% if ( save_flag == 1 )
%     fid = fopen(f,'w');
%     [nrows,ncols] = size(table);
%     
%     for row = 1:nrows
%         fprintf(fid,'%s\n',table{row,:});
%     end
%     fclose(fid);
% end

% % figure(1)
% % print('-depsc2', ['Helm1D_',bc,'_FS_timing_in_k.eps'])
% % figure(2)
% % print('-depsc2', ['Helm1D_',bc,'_FS_iterations_in_k.eps'])
%
% break

% % Plot the relative residual curves FOR ONE WAVENUMBER k
% figure(1)
% semilogy(0:iter_lu(2), resvec_lu/resvec_lu(1), 'k-')
% hold on
% semilogy(0:iter_ilu(2), resvec_ilu/resvec_ilu(1), 'k-')
% semilogy(0:iter_mg(2), resvec_mg/resvec_mg(1), 'k-')
% semilogy(0:iter_ex_FS(2), resvec_ex_FS/resvec_ex_FS(1), 'b-')
% semilogy(0:iter_lu_FS(2), resvec_lu_FS/resvec_lu_FS(1), 'b-')
% semilogy(0:iter_ilu_FS(2), resvec_ilu_FS/resvec_ilu_FS(1), 'b-')
% semilogy(0:iter_mg_FS(2), resvec_mg_FS/resvec_mg_FS(1), 'b-')
% hold off
% ylabel('relative residual')
% xlabel('iteration')
% legend('f_{ex} = M^{-1}A*x (plain gmres)', 'f_{lu}', 'f_{ilu}', ...
%     'f_{ex} and poly acc with FS', 'f_{lu} and poly acc with FS', ...
%     'f_{ilu} and poly acc with FS')
% title(['1D Helmholtz with CSL-preconditioner (k=',num2str(k),')']) FS = 22 ; % font size
% set(gca,'LooseInset',get(gca,'TightInset'))
% set(gca,'FontSize',FS);
%
% print('-depsc2', 'Helm1D_rr_FS_different_MinvA.eps')
% table_data = [wavenum', ...,
%               time_lu, iter_lu(:,2),...
%               time_lu_FS, iter_lu_FS(:,2),...
%               time_mg, iter_mg(:,2),...
%               time_mg_FS, iter_mg_FS(:,2)];

% name_table = strcat('table_time_vs_k_',dimn,'D_',...
%                  bdc,'_','kmin',kmin,'_',...
%                  'kmax',kmax,'.tex');
% clear input;
% input.data = table_data;
%
% %We want a complete LaTex document (change this later)
% input.makeCompleteLatexDocument = 1;
% input.dataFormatMode = 'column'; % use 'column' or 'row'. if not set 'colum' is used
%
% %Set column labels (use empty string for no label):
% input.tableColLabels = {'k','Time Exact CSL (s)','Iter. Exact CSL (s)',...
%                         'Time Faber with Exact CSL (s)','Iter. Faber series with Exact CSL',...
%                         'Time MG CSL (s)','Iter. MG-CSL',...
%                         'Time Faber series with MG-CSL(s)','Iter. Faber series with MG CSL'};
%
% %Set row labels (use empty string for no label):
% %input.tableRowLabels = {'','','','','','','','','',''};
%
% % Column alignment in Latex table ('l'=left-justified, 'c'=centered,'r'=right-justified):
% input.tableColumnAlignment = 'c';
%
% % Switch table borders on/off (borders are enabled by default):
% input.tableBorders = 1;
%
% % Formatting-string to set the precision of the table values:
% % For using different formats in different rows use a cell array like
% % {myFormatString1,numberOfValues1,myFormatString2,numberOfValues2, ... }
% % where myFormatString_ are formatting-strings and numberOfValues_ are the
% % number of table columns or rows that the preceding formatting-string applies.
% % Please make sure the sum of numberOfValues_ matches the number of columns or
% % rows in input.tableData!
% %
% input.dataFormat = {'%.0f',1,'%.2f',1,'%.0f',1,'%.2f',1,'%.0f',1,'%.2f',1,'%.0f',1,'%.4f',1,'%.0f',1,}; % three digits precision for first two columns, one digit for the last
%
% % Switch table borders on/off (borders are enabled by default):
% input.tableBorders = 1;
%
% %Switch to generate a complete LaTex document or just a table:
% input.makeCompleteLatexDocument = 1;
%
% %Generate LaTex code
% latex = latexTable(input);
%
% % save LaTex code as file
% fid=fopen(name_table,'w');
% [nrows,ncols] = size(latex);
% for row = 1:nrows
%     fprintf(fid,'%s\n',latex{row,:});
% end
% fclose(fid);
%
% %fprintf('\n... your LaTex code has been saved as ''Table1.tex'' in your working directory\n');

