%Field of values of finite element matrices 1D
%ADEF preconditioner
%
% Experiments with the  Helmholtz equation.
% preconditoned with the shifted Laplacian (SL)
% and the deflated SL
%
% We compare iteration numbers for solving with GMRES
% AAeps^{-1}x = b and
% ABx=b, where B is the deflated shifted Laplacian
% 
%
%% Construction of the matrices
close all

save_flag = 1;  % save_flag=1: save plots and table, =0 do not save.

% Setup parameters
% Setup list of wavenumbers and shifts
dim = 2;
poweps    = 2;
factoreps = [0.5 1];

%Wavenumber
kk      = [10 20 40 60];
%kk = 10;

bc = 'som';
plot_defcsl = 'yes';
plot_csl_defcsl = 'yes';

%Parameters for GMRES
restart   = [];
tol       = 1e-8;
maxit     = 200;

linetyp = {':','-.','--'};
color   = [0 0 0; 0.5 0 0.5; 0 0 1; 0 0.5 0; 1 0 0];

ppw = 0.5;
%if pollution = 'no' the number of points np= ceil(k^(3/2))
pollution = 'no';
iter_adef = zeros(length(kk),length(factoreps));
iter_csl  = zeros(length(kk),length(factoreps));

%% Comparison of GMRES with Shifted Laplace and deflated SL problems
for i=1:length(kk)
    k   = kk(i);  
    for j=1:length(factoreps)
        eps = factoreps(j)*k^poweps;
        
        %choosing the number of points
        npf = ceil(ppw*k/(2*pi))-1;
        if strcmp(pollution,'no')
            npf = ceil(k^(3/2));
        end
        
        if (mod(npf+1,2)==1)  %set an even number of interior points in 1D
            npf = npf+1;
        end
        
        npc  = (npf-1)/2;
        H = 1/(npc+1);
          
        [node1,elem1] = squaremesh([0,1,0,1],H);  %coarse mesh
        [node,elem] = uniformrefine(node1,elem1); %fine mesh (uniformly refined)
        
        %Find boundary nodes
        [bdNode,bdEdge,isBdNode] = findboundary(elem);
        [bdNode1,bdEdge1,isBdNode1] = findboundary(elem1);

        
        %Sets Sommerfeld boundary conditions on all boundary edges
        bdFlag = setboundary(node,elem,'ABC');
        bdFlag1 = setboundary(node1,elem1,'ABC');

        %The structures pde(helm,SL) contain data for the Helmholtz and
        %shifted Laplace problems
        pdehelm = helmholtz2Dconstantwndata(k,0,1);
        pdeSL   = helmholtz2Dconstantwndata(k,factoreps(j),poweps);
        
        option.tol = 1e-8;
        fprintf('beginning computation of fem matrices for k=%d  \n', k);
        tic
        [eqn1,~] = helmholtz2Dfem(node,elem,pdehelm,bdFlag,bdEdge);
        [eqn2,~] = helmholtz2Dfem(node,elem,pdeSL,bdFlag,bdEdge);
        [eqncoarse,~] = helmholtz2Dfem(node1,elem1,pdehelm,bdFlag1,bdEdge1);
    
  
        %Helmholtz and shifted Laplace matrices
        A    = eqn1.A;
        Aeps = eqn2.A;
        Acoarse = eqncoarse.A; 
        
        
        fprintf('finished computation of fem matrices  for k=%d  \n', k);
        fprintf('size of fem matrices for k=%d: %d \n', k, length(eqn1.A));
        
  
        option.twolevel = true;
        numlev = 2;
        [mg_mat,~,R,~] = mg_setupfem_2D(npc,numlev,pdehelm,option);
          
        
        Z = R'; 
        fprintf('beginning computation of lu of Aeps for k=%d  \n', k);
        [L,U]    = lu(Aeps);
        
        
        Ac       = Z'*A*Z; %Coarse operator
        fprintf('beginning computation of lu of Ac for k=%d  \n', k);
        
        [Lc, Uc] = lu(Ac);
        
        %Inverses of SLaplace and coarse operator
        N        = length(A);
        I        = speye(N);
        Aepsinv  =  @(x) U\(L\x);     %Inverse of shifted Laplace
        Acinv    =  @(x) Uc\(Lc\x);   %Inverse of coarse Helmholtz
        Q        =  @(x) Z*Acinv(Z'*x);
        
        %Definition of ADEF preconditioner
        B = @(x) Aepsinv(x-A*Q(x))+Q(x);
        
        %Definition of preconditioned operators
        AAepsinv  = @(x) A*(U\(L\x));  %SL
        AB =  @(x) A*B(x);  %Deflated SL
        
        %Input for GMRES
        b  = ones(length(A),1);
        x0 = zeros(size(b));
        
        %GMRES runs
        fprintf('beginning adef gmres run for k=%d  \n', k);
        
        [~, ~, ~, iter, ~] = gmres(AB, b, restart, tol, maxit);
        iter_adef(i,j) = iter(2);
        
        if strcmp(plot_csl_defcsl,'yes');
            fprintf('beginning csl gmres run for k=%d  \n', k);
            [~, ~, ~, iter, ~] = gmres(AAepsinv, b, restart, tol, maxit);
            iter_csl(i,j) = iter(2);
        end
        
    end %End of factoreps for-loop
end %end of wavenum loop


%% Plot iteration numbers and save as tikz (.tex) and .eps
figure(1)

%plot iteration counts
%preconditioned problem
plot(kk, iter_adef(:,1),'*','Color',color(1,:),'Markersize',10);
hold on
if strcmp(plot_csl_defcsl, 'yes')
    plot(kk, iter_csl(:,1),'Color',color(1,:),...
        'linestyle','-','Linewidth',3);
end

%preconditioned problem
for j=2:length(factoreps)
    plot(kk, iter_adef(:,j),'*','Color',color(j,:));
    
    if strcmp(plot_csl_defcsl,'yes')
        plot(kk, iter_csl(:,j),'Color',color(j,:),...
            'linestyle','-','Linewidth',3);
    end
    
end

axis([10 200 0 110]);

%set(gca,'Xtick',[0 0.5 1],'FontSize',30);
%set(gca,'Ytick',[-0.5 0 0.5],'FontSize',30);

hold off
ylabel('Number of GMRES iterations','fontsize',25)
xlabel('Wavenumber','fontsize',25)
%legend('No Prec.', ...
legend('\epsilon=k^2', ...
    '\epsilon=2k^2',...
    '\epsilon=10k^2',...
    'Location','NorthEast')
%title(['1D Helmholtz with CSL-preconditioner and poly. accel., deg=',...
%num2str(deg)])
FS = 25; % font size
set(gca,'LooseInset',get(gca,'TightInset'))
set(gca,'FontSize',FS)

plot_iter_tex = strcat('comp_iter_vs_k_eps_1D_som.tex');

%get path of current .m file
currentpath  = mfilename('fullpath');
%generate file name, for saving in tex_files/figures/new_exp
plot_file_tex = fullfile(currentpath,plot_iter_tex);

if ( save_flag == 1 )
    %Save as tikz figure in .tex file
    matlab2tikz('filename',plot_iter_tex,'standalone',true,'extraaxisoptions',...
        ['xlabel style={font=\LARGE},', ...
        'ylabel style={font=\LARGE},']);
    
    %Save the plot as .eps
    % plot_iter_eps = strcat('comp_iter_vs_k_eps_1D_som.eps');
    % plot_file_eps = fullfile(path,plot_iter_eps);
    % print('-depsc2',plot_file_eps)
end
