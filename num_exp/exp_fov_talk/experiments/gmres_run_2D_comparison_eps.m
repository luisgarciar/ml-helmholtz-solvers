%% Experiments with the  Helmholtz equation.
% preconditoned with the shifted Laplacian (SL)
% and the deflated SL
% We compare iteration numbers for solving
% A inv(Aeps)x = b
%
%% Construction of the matrices
clear all
close all
save_flag = 1;  % save_flag=1: save plots and table, =0 do not save.

% Setup parameters
% Setup list of wavenumbers and shifts
dim       = 2;
poweps    = [0 1 1.5 2];
%poweps   = [0];
%poweps   = [1.5];



factoreps = 1;

%Wavenumber
k = 40;

bc = 'som';
plot_defcsl = 'yes';
plot_csl_defcsl = 'yes';

%Parameters for GMRES
restart   = [];
tol       = 1e-8;
maxit     = 150;

linetyp = {':','-.','--'};
color   = [0 0 0; 0.5 0 0.5; 0 0 1; 0 0.5 0; 1 0 0];

ppw = 0.5;
%if pollution = 'no' the number of points np= ceil(k^(3/2))
pollution = 'no';
iter_csl  = zeros(1,length(factoreps));


%% GMRES comparison of shifted Laplace problems
for j=1:length(poweps)
    eps = factoreps*k^poweps(j);
    
    %choosing the number of points
    if strcmp(pollution,'no')
        npf = ceil(k^(3/2));
    end
    
    if (mod(npf,2)==0)  %set an odd number of interior points in 1D
        npf = npf+1;
    end
    npc = (npf-1)/2;
    
    bc = 'som';
    
    %npf : number of interior points
    %Construct coarse and fine grid meshes of meshsizes H,h
    H = 1/(npc+1);
    
    [node,elem] = squaremesh([0,1,0,1],H);  %coarse mesh
    [node,elem] = uniformrefine(node,elem); %fine mesh (uniformly refined)
    
    %Find boundary nodes
    [bdNode,bdEdge,isBdNode] = findboundary(elem);
    
    %Sets Sommerfeld boundary conditions on all boundary edges
    bdFlag = setboundary(node,elem,'ABC');
    
    %The structures pde(helm,SL) contain data for the Helmholtz and
    %shifted Laplace problems
    pdehelm    = helmholtz2Dconstantwndata(k,0,1);
    pdeSL      = helmholtz2Dconstantwndata(k,factoreps,poweps(j));
    option.tol = 1e-12;
    
    fprintf('beginning computation of fem matrices for k=%d  \n', k);
    tic
    [eqn1,~] = helmholtz2Dfem(node,elem,pdehelm,bdFlag,bdEdge);
    [eqn2,~] = helmholtz2Dfem(node,elem,pdeSL,bdFlag,bdEdge);
    time_fem = toc;
    
    fprintf('finished computation of fem matrices  for k=%d  \n', k);
    fprintf('time fem: %f \n', time_fem);
    fprintf('size of fem matrices for k=%d: %d \n', k, length(eqn1.A));
    
    %Helmholtz and shifted Laplace matrices
    A    = eqn1.A;
    Aeps = eqn2.A;
    
    %Inverses of SLaplace and coarse operator
    [L,U]   = lu(Aeps);
    Aepsinv =  @(x) U\(L\x);     %Inverse of shifted Laplace
    
    %Definition of preconditioned operators
    AAepsinv = @(x) A*Aepsinv(x);  %SL
    
    %Input for GMRES
    b  = ones(length(A),1);
    x0 = zeros(size(b));
    
    %GMRES runs
    if isequal(poweps(j),0)
        [~, flag, relres, iter,resvec] = gmres(A, b, restart, tol, maxit);
    else
        [~, flag,relres, iter, resvec] = gmres(AAepsinv, b, restart, tol, maxit);
    end
    
    semilogy(1:(iter(2)+1),resvec/norm(b),'-','Color',color(j,:),'Markersize',10,'Linewidth',2);
    xlabel('Iteration number','fontsize',20)
    ylabel('Relative Residual Norm','fontsize',20)
    FS = 25; % font size
    set(gca,'LooseInset',get(gca,'TightInset'))
    set(gca,'FontSize',FS)
   
    hold on
end %End of factoreps for-loop

leg1 = legend('No prec.','$\varepsilon=k$','$\varepsilon=k^{1.5}$','$\varepsilon=k^2$');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',20);

% %% Plot iteration numbers and save as tikz (.tex) and .eps
% figure(1)
%
% %plot iteration counts
% %preconditioned problem
% plot(kk, iter_adef(:,1),'*','Color',color(1,:),'Markersize',10);
% hold on
%
% %preconditioned problem
% for j=2:length(factoreps)
%     plot(kk, iter_adef(:,j),'*','Color',color(j,:),'Markersize',10);
% end
%
% axis([20 60 0 10]);
%
% %set(gca,'Xtick',[0 0.5 1],'FontSize',30);
% %set(gca,'Ytick',[-0.5 0 0.5],'FontSize',30);
%
% hold off
% ylabel('Number of GMRES iterations','fontsize',20)
% xlabel('Wavenumber','fontsize',20)
% %legend('No Prec.', ...
% legend('\epsilon=k^2', ...
%     '\epsilon=2k^2',...
%     '\epsilon=10k^2',...
%     'Location','NorthEast')
% %title(['1D Helmholtz with CSL-preconditioner and poly. accel., deg=',...
% %num2str(deg)])
% FS = 25; % font size
% set(gca,'LooseInset',get(gca,'TightInset'))
% set(gca,'FontSize',FS)
%
% plot_iter_tex = strcat('comp_iter_vs_k_eps_1D_som.tex');
%
% %get path of current .m file
% currentpath  = mfilename('fullpath');
% %generate file name, for saving in tex_files/figures/new_exp
% plot_file_tex = fullfile(currentpath,plot_iter_tex);
%
% if ( save_flag == 1 )
%     %Save as tikz figure in .tex file
%     matlab2tikz('filename',plot_iter_tex,'standalone',true,'extraaxisoptions',...
%         ['xlabel style={font=\LARGE},', ...
%         'ylabel style={font=\LARGE},']);
%
%     %Save the plot as .eps
%     % plot_iter_eps = strcat('comp_iter_vs_k_eps_1D_som.eps');
%     % plot_file_eps = fullfile(path,plot_iter_eps);
%     % print('-depsc2',plot_file_eps)
% end
