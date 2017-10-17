%Field of values of finite element matrices 1D
%ADEF preconditioner

%% Construction of the matrices
clear all
close all

dim = 2;
poweps    = 2;
factoreps = 1;
bc = 'som';

%Wavenumber
%kk = 200;
kk = 50;
%kk      = [20 50 100 150];
iter_SL = zeros(length(kk),2);

%Colors for plots
hexcolor   = ['#332288'; '#88CCEE'; '#44AA99'; '#117733'; '#999933'; ...
    '#DDCC77'; '#CC6677'; '#882255'; '#AA4499'];
rgbcolor1  = hex2rgb(hexcolor);
%rgbcolor2 = linspecer(9);

minfov  = zeros(length(kk),1);
linetyp = {'-','.','--','-'};
color   = {'m','b','g','k'};
marker  = {'*','o','.','x'};
opt     = {'m','b','g','k'};
%str2={'k','k.-','k--','b*-'};

plot_csl     = 'yes';
plot_gmres   = 'no';
ppw = 0.5;

%if pollution = 'no' the number of points np= ceil(k^(3/2))
pollution = 'no';
fvpts = 50;

%gmres parameters
restart = [];
tol     = 1e-8;
maxit   = 100;


%% Plot of FOV of Shifted Laplace problems
for i=1:length(kk)
    k   = kk(i);
    eps = factoreps*k^poweps;
    
    %choosing the number of points
    npf = ceil(ppw*k/(2*pi))-1;
    
    if strcmp(pollution,'no')
        npf = ceil(k^(3/2));
    end
    
    if (mod(npf,2)==0)  %set an odd number of interior points in 1D
        npf = npf+1;
    end
    npc = (npf-1)/2;
    
    H = 1/(npc+1);
                              
    [node,elem] = squaremesh([0,1,0,1],H);  %coarse mesh
    [node,elem] = uniformrefine(node,elem); %fine mesh (uniformly refined)
    
    %Find boundary nodes
    [bdNode,bdEdge,isBdNode] = findboundary(elem);
    
    %Sets Sommerfeld boundary conditions on all boundary edges
    bdFlag = setboundary(node,elem,'ABC');
    
    %The structures pde(helm,SL) contain data for the Helmholtz and
    %shifted Laplace problems
    pdehelm = helmholtz2Dconstantwndata(k,0,1);
    pdeSL   = helmholtz2Dconstantwndata(k,factoreps,poweps);
    
    option.tol = 1e-12;
    fprintf('beginning computation of fem matrices for k=%d  \n', k);
    tic
    [eqn1,~] = helmholtz2Dfem(node,elem,pdehelm,bdFlag,bdEdge);
    [eqn2,~] = helmholtz2Dfem(node,elem,pdeSL,bdFlag,bdEdge);
    time_fem = toc;

    fprintf('finished computation of fem matrices  for k=%d  \n', k);
    fprintf('time fem: %f  \n', time_fem);
    fprintf('size of fem matrices for k=%d: %d \n', k, length(eqn1.A));
    
    %Helmholtz and shifted Laplace matrices
    A    = eqn1.A;
    Aeps = eqn2.A;
    
    %We compute FOV(MAepsinv), since FOV(AAepsinv) = 1+ FOV(MAepsinv)
    %fprintf('press any key to continue \n');
    %pause;
    tic
    fprintf('beginning computation of lu of Aeps for k=%d  \n', k);
    [L,U] = lu(Aeps);
    time_lu = toc;
    fprintf('lu factorization for k=%d finished \n', k);
    fprintf('time lu: %f  \n', time_lu);

    LH = U'; UH = L';
    N  = length(A);
    
    %Let C = MAepsinv
    M = eqn1.M;
    C    = @(x) M*(U\(L\x));
    CH   = @(x) UH\(LH\(M*x));
    HerC = @(x) 0.5*(feval(C,x) + feval(CH,x)); %Hermitian part of C
    
    %Compute first the maximum eigenvalue of A
%   [eigvecA, eigvalA] =  eig(full(A));
%   [eigvalmaxA,ind] = max(diag(eigvalA));
%   vmaxA = eigvecA(:,ind);
         
    %Compute the maximum eigenvalue of the Hermitian part of A
    opts.isreal = 0;
    opts.p      = 30;
    opts.issym  = true;
    opts.tol    = 1e-6;  
    opts.v0     = rand(N,1);
    fprintf('beginning computation of max eigvalue of H(C) for k = %d  \n', k);
    tic
    [vmaxHerC, eigmaxHerC] =  eigs(HerC,N,1,'LM',opts);
    time_maxeigv = toc;
    fprintf('computation of max eigenvalue of H(C) for k=%d finished \n', k);
    fprintf('time maxeigv: %f  \n', time_maxeigv);
    
    %Use the maximum eigenvalue of the Hermitian part of A
    %to initiate computation of fovC
    fprintf('beginning computation of fov for k=%d \n', k);
    tic
    [fovC,~,~] = sfov(C,CH,vmaxHerC,N,50);
    time_fov = toc;
    fprintf('time fov: %f  \n', time_fov);
    
    fovAhat = 1+1i*eps*fovC;
    reFOV  = real(fovAhat); imFOV = imag(fovAhat);
    
    cvh = convhull(reFOV,imFOV);
    
    %plotting the fov
    label       = strcat('k = ',num2str(k));
    fovplot(i)  = plot(reFOV(cvh),imFOV(cvh),'Color',rgbcolor1(i,:),...
        'LineWidth',3,'DisplayName',label);
    hold on
    plot(0,0,'k*','Markersize',10,'LineWidth',2);
    axis equal
    axis([-0.2 1.2 -0.7 0.7]);
    xlabel('Re(z)','FontSize',18);
    ylabel('Im(z)','FontSize',18);
    set(gca,'Xtick',[0 0.5 1],'FontSize',18);
    set(gca,'Ytick',[-0.5 0 0.5 1],'FontSize',18);
    
end


legend(fovplot);

%Filename format:
%1d_fov_csl_kmin_kmax_pointswavelength_realshift_imagshift.tex

kmin  = num2str(min(kk));
kmax  = num2str(max(kk));
pts   = num2str(ppw);
powershift  = num2str(poweps);
factorshift = num2str(10*factoreps);


%Tikz Axis formatting
x = xlabel('$\mathrm{Re}(z)$');
% x-axis label
set(x,'Interpreter','latex')
y=ylabel('$\mathrm{Im}(z)$','interpreter','latex'); % x-axis label
set(y,'Interpreter','latex')
figure(1)
name1 = strcat('1d_fov_csl_kmin',kmin,'_kmax',kmax,'_ppw',pts, ...
    '_pshift_',powershift,'_fshift_',factorshift,'.tex');

if strcmp(pollution,'no')
    name1 = strcat('1d_fov_csl_kmin',kmin,'_kmax',kmax,'nopoll',pts, ...
        '_pshift_',powershift,'_fshift_',factorshift,'.tex');
end
matlab2tikz('filename',name1,'standalone',true); %save .tex file of tikz figure%
%close all;


if(strcmp(plot_gmres,'yes'))
    
    for i=1:length(kk)
        k   = kk(i);
        eps = factoreps*k^poweps;
        
        %choosing the number of points
        npf = ceil(ppw*k/(2*pi))-1;
        if strcmp(pollution,'no')
            npf = ceil(k^(3/2));
        end
        
        if (mod(npf+1,2)==1)  %set an even number of interior points in 1D
            npf = npf+1;
        end
        npc = (npf-1)/2;
        
        A    = helmholtzfem(k,npf,0,bc);           %Helmholtz matrix
        Aeps = helmholtzfem(k,npf,eps,bc);         %Shifted Laplace matrix
        b = ones(length(A),1);
        
        [L, U] = lu(Aeps);
        Ahat   = @(x) A*(U\(L\x));
        [~,~,~,iter_SL(i,:),resvec] = gmres(Ahat,b,restart,tol,maxit);
        
    end
    
    close all
    figure(2)
    plot(kk, iter_SL(:,2), 'k+')
    hold on
    
end



