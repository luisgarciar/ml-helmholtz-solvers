%Field of values of finite element matrices 1D
%ADEF preconditioner

%% Construction of the matrices
close all

dim       = 1;
poweps    = 2;
factoreps = 5;
bc = 'som';

%Wavenumber
kk  = 20;
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
ppw          =  0.5;

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
    
    %% Sparse FOV of Helmholtz + Shifted Laplace
    tic
    fprintf('beginning computation of lu of Aeps for k=%d  \n', k);
    [L,U] = lu(Aeps);
    UH = L'; LH=U';
    time_lu = toc;
    fprintf('lu factorization for k=%d finished \n', k);
    fprintf('time lu: %f  \n', time_lu);
    
    N  = length(A);
    
    Ahat    = @(x) A*(U\(L\x));
    AhatH   = @(x) UH\(LH\(A'*x));
    H       = @(x) 0.5*(feval(Ahat,x) + feval(AhatH,x)); %Hermitian part of Ahat
    
    
    %% Sparse FOV of Deflated shifted Laplacian
    
    %Set up two-level structure
    option.twolevel = true;
    numlev = 2;
    [mg_mat,~,~,Z] = mg_setupfem_2D(npc,numlev,pdehelm,option);
    
    %Deflation subspace: columns of Z
    dim_def = size(Z,2);
    
    %Deflated-shifted Operator PadefA
    M  = eqn1.M;
    Ac = Z'*A*Z; %Coarse operator
    
    %We compute the projection P = I-AZ(Z'AZ)^(-1)Z
    %using the procedure from Stewart,
    %"On the Numerical Analysis of Oblique Projectors" - SIMAX, 2011
    
    fprintf('Beginning computation of qr for projections for k=%d  \n', k);
    tic
    [X,S] = qr(full(A*Z));
    time_qr1 = toc;
    fprintf('finished computation of qr1 for k=%d  \n', k);
    fprintf('time qr1: %f  \n', time_qr1);
 
    tic
    [Y,T] = qr(full(Z));
    time_qr2 = toc;
    fprintf('finished computation of qr2 for k=%d  \n', k);
    fprintf('time qr2: %f  \n', time_qr2);
    
    [L1,U1] = lu(Y'*X);
    [L2,U2] = lu(X'*Y);
        
%   M1 = Y'*X; 
%   tic; [Q1,R1] = qr(M1); time_qr3 = toc;
%   fprintf('finished computation of qr3 for k=%d  \n', k);
%   fprintf('time qr3: %f  \n', time_qr3);
     
%   M2 = X'*Y; 
%   tic; [Q2,R2] = qr(M2); time_qr4 = toc;
%   fprintf('finished computation of qr for k=%d  \n', k);
%   fprintf('time qr3: %f  \n', time_qr4);
 
    fprintf('finished computation of qr for k=%d \n', k);
    %fprintf('total time qr: %f \n', time_qr1+time_qr2+time_qr3+time_qr4);
      
% We compute the projection  P = AZ(Z'AZ)^(-1)AZ

%   P  = @(x) X*(R1\(Q1'*(Y'*x)));
%   PH = @(x) Y*(R2\(Q2'*(X'*x)));     
        
    P  = @(x) X*(U1\(L1\(Y'*x)));
    PH = @(x) Y*(U2\(L2\(X'*x)));

% fprintf('beginning computation of lu of coarse matrix for k=%d  \n', k);
% tic
% [Lc, Uc] = lu(Ac);
% time_luc = toc;

% fprintf('coarse lu factorization for k=%d finished \n', k);
% fprintf('time lu coarse: %f  \n', time_luc);
%
% LcH = Uc'; UcH = Lc';
% sqrtM = sqrtm(full(M));

    
    %Deflated operator
    I        =  speye(N);
    Aepsinv  =  @(x) U\(L\x);     %Inverse of shifted Laplace
    AepsHinv =  @(x) UH\(LH\x);   %Inverse of Hermitian transpose of shifted Laplace
       
    %Acinv    =  @(x) Uc\(Lc\x);   %Inverse of coarse Helmholtz
    %AcHinv   =  @(x) UcH\(LcH\x); %Inverse of Hermitian transpose of coarse Helmholtz
    
    %FOV in the Euclidean inner product
    %AP  =  @(x)  M*Aepsinv(x-A*Z*Acinv((Z'*x)));
    %APH =  @(x)  AepsHinv(M*x)- Z*AcHinv(Z'*A'*AepsHinv(M*x));
        
    %Field of values of APadef in Minv inner product
    %AP  =    @(x)  sqrtM*Aepsinv((sqrtM*x-A*Z*Acinv((Z'*sqrtM*x))));
    %APH =    @(x)  sqrtM*(AepsHinv(sqrtM*x)- Z*AcHinv(Z'*A'*AepsHinv(sqrtM*x)));
        
    %FOV using the XQRY form of the projection
    AP  = @(x) M*Aepsinv(x-feval(P,x));
    APH = @(x) AepsHinv(M*x) - feval(PH,AepsHinv(M*x));
    
    fprintf('beginning computation of fov for k=%d \n', k);
    tic
    vmaxH  = rand(N,1);
    [fovAP1,~,~] = parallel_sfov(AP,APH,vmaxH,N,50);
    fovAP1       = 1+1i*eps*fovAP1;
    time_fov = toc;
    fprintf('\n time fov: %f  \n', time_fov);
    
    %close all
    refovAP1 = real(fovAP1); imfovAP1= imag(fovAP1); 
    figure(1)
    plot(refovAP1,imfovAP1,'b','LineWidth',3)  %Plot the field of values
    hold on
    plot(0,0,'k+','Markersize',10,'LineWidth',2);
    axis('equal');
    axis([-0.2 1.2 -0.7 0.7]);
    xlabel('Re(z)','FontSize',14);
    ylabel('Im(z)','FontSize',14);
    set(gca,'Xtick',[ 0   0.5   1],'FontSize',14);
    set(gca,'Ytick',[-0.5 0 0.5 1],'FontSize',14);
    
    
    x = xlabel('$\mathrm{Re}(z)$');     % x-axis label
    set(x,'Interpreter','latex');
    y=ylabel('$\mathrm{Im}(z)$','interpreter','latex'); % x-axis label
    set(y,'Interpreter','latex');
    figure(1);
    
    wn  = num2str(kk(i));
    pts   = num2str(ppw);
    powershift  = num2str(poweps);
    factorshift = num2str(10*factoreps);
    
    
    name1 = strcat('2D_fov_dcsl_wn',wn,'_ppw',pts, ...
        '_pshift_',powershift,'_fshift_',factorshift,'.tex');
    
    if strcmp(pollution,'no')
        name1 = strcat('2D_fov_dcsl_wn',wn,'_nopoll', ...
            '_pshift_',powershift,'_fshift_',factorshift,'.tex');
    end
    
    matlab2tikz('filename',name1,'standalone',true); %save .tex file of tikz figure%
    %End of FOV DCSL 2D
    
    
end
