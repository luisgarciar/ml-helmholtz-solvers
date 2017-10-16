%Field of values of finite element matrices 1D
%ADEF preconditioner

%% Construction of the matrices
clear all
close all

dim = 1;
poweps    = 2;
factoreps = 1;
bc = 'som';

%Wavenumber
kk  = 10;
%kk = [20 50 100 150];
iter_SL = zeros(length(kk),2);

%Colors for plots
hexcolor   = ['#332288'; '#88CCEE'; '#44AA99'; '#117733'; '#999933'; ...
              '#DDCC77'; '#CC6677'; '#882255'; '#AA4499'];
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
fvpts = 50;

%gmres parameters
restart = [];
tol     = 1e-8;
maxit   = 100;

ppw = 0.5;
%if pollution = 'no' the number of points np= ceil(k^(3/2))
pollution = 'no';

for i=1:length(kk)
    close all
    
    k   = kk(i);
    eps = factoreps*k^poweps;
    npf = ceil(k^(3/2));
    
    if (mod(npf,2)==0)  %set an odd number of interior points in 1D
        npf = npf+1;
    end
    npc = (npf-1)/2;
    
    poweps    = 1;
    factoreps = 1;
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
    pdehelm = helmholtz2Dconstantwndata(k,0,1);
    pdeSL   = helmholtz2Dconstantwndata(k,factoreps,poweps);
    
    option.tol = 1e-12;
    [eqn1,~] = helmholtz2Dfem(node,elem,pdehelm,bdFlag,bdEdge);
    [eqn2,~] = helmholtz2Dfem(node,elem,pdeSL,bdFlag,bdEdge);
    
    %Helmholtz and shifted Laplace matrices
    A    = eqn1.A;
    Aeps = eqn2.A;
    
    %% Sparse FOV of Helmholtz + Shifted Laplace
    [L,U] = lu(Aeps);
    LH = U'; UH = L';
    N  = length(A);
    
    Ahat    = @(x) A*(U\(L\x));
    AhatH   = @(x) UH\(LH\(A'*x));
    H       = @(x) 0.5*(feval(Ahat,x) + feval(AhatH,x)); %Hermitian part of Ahat
    
    %Compute first the maximum eigenvalue of A
    [eigvecA, eigvalA] =  eig(full(A));
    [eigvalmaxA,ind] = max(diag(eigvalA));
    vmaxA = eigvecA(:,ind);
    
    %Compute now the maximum eigenvalue of the Hermitian part of A
    opts.p = 30;
    opts.isreal = 0;
    opts.v0 = vmaxA;
    [vmaxH, eigmaxH] =  eigs(H,N,1,'LM',opts);
    
    %Use the maximum eigenvalue of the Hermitian part of A
    [fovAhat,~,~] = sfov(Ahat,AhatH,vmaxH,N,50);
    
    %Find the convex hull of the FOV and plot it
     reFOV = real(fovAhat); imFOV = imag(fovAhat);
      cvh = convhull(reFOV,imFOV);
    
    %% Sparse FOV of Deflated shifted Laplacian
    
    %Set up two-level structure 
    option.twolevel = true;
    numlev = 2;
    [mg_mat,~,R,Z] = mg_setupfem_2D(npc,numlev,pdehelm,option);
    
   %Deflation subspace: columns of Z
    dim_def = size(Z,2);
 
    %Deflated-shifted Operator PadefA
    M = eqn1.M;
    Ac = Z'*A*Z; %Coarse operator
    [Lc, Uc] = lu(Ac);
    LcH = Uc'; UcH = Lc';
    %sqrtM = sqrtm(full(M));
    
    %Deflated operator
    I        = speye(N);
    Aepsinv  =  @(x) U\(L\x);     %Inverse of shifted Laplace
    AepsHinv =  @(x) UH\(LH\x);   %Inverse of Hermitian transpose of shifted Laplace
    Acinv    =  @(x) Uc\(Lc\x);   %Inverse of coarse Helmholtz
    AcHinv   =  @(x) UcH\(LcH\x); %Inverse of Hermitian transpose of coarse Helmholtz
    
   %Field of values of APadef in Minv inner product
    %AP  =    @(x)  sqrtM*Aepsinv((sqrtM*x-A*Z*Acinv((Z'*sqrtM*x))));
    %APH =    @(x)  sqrtM*(AepsHinv(sqrtM*x)- Z*AcHinv(Z'*A'*AepsHinv(sqrtM*x)));
    
    %FOV in the Euclidean inner product
    AP  =  @(x)  M*Aepsinv(x-A*Z*Acinv((Z'*x)));
    APH =  @(x)  AepsHinv(M*x)- Z*AcHinv(Z'*A'*AepsHinv(M*x));
  
    [fovAP,~,~] = sfov(AP,APH,vmaxH,N,50);
    fovAP       = 1+1i*eps*fovAP;    
    %close all
   
    refovAP = real(fovAP); imfovAP= imag(fovAP);
    plot(refovAP,imfovAP,'k','LineWidth',1)      % Plot the field of values
    hold on
    plot(0,0,'k+','Markersize',10,'LineWidth',2)
    axis('equal');
    axis([-0.2 1.2 -0.7 0.7]);
    xlabel('Re(z)','FontSize',14);
    ylabel('Im(z)','FontSize',14);
    set(gca,'Xtick',[ 0   0.5   1],'FontSize',14);
    set(gca,'Ytick',[-0.5 0 0.5 1],'FontSize',14);
    
    x = xlabel('$\mathrm{Re}(z)$');
    % x-axis label
    set(x,'Interpreter','latex')
    y=ylabel('$\mathrm{Im}(z)$','interpreter','latex'); % x-axis label
    set(y,'Interpreter','latex')
    figure(1)
    
    wn  = num2str(kk(i));
    pts   = num2str(ppw);
    powershift  = num2str(poweps);
    factorshift = num2str(10*factoreps);
 
    
    name1 = strcat('1d_fov_dcsl_wn',wn,'_ppw',pts, ...
        '_pshift_',powershift,'_fshift_',factorshift,'.tex');
    
    if strcmp(pollution,'no')
        name1 = strcat('1d_fov_dcsl_wn',wn,'_nopoll', ...
            '_pshift_',powershift,'_fshift_',factorshift,'.tex');
    end
    
    matlab2tikz('filename',name1,'standalone',true); %save .tex file of tikz figure%
    %End of FOV DCSL
    
    
end
