%Field of values of finite element matrices 1D
%ADEF preconditioner

%% Construction of the matrices
close all

dim       = 1;
poweps    = 2;
factoreps = 5;
bc = 'som';

%Wavenumber
kk  = 40;
kk = [10 20 40];
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
    
    %npf : number of interior points in fine grid
    npf = ceil(k^(3/2));
    if (mod(npf,2)==0)  %set an odd number of interior points in 1D
        npf = npf+1;
    end
    npc = (npf-1)/2; %number of interior points in coarse grid
    
    poweps    = 1;
    factoreps = 1;
    bc        = 'som';
    dim       = 2;
    
    %Construct coarse and fine grid meshes of meshsizes h, H=2h
    numlev  = 2;
    op_type = 'gal';
    [mg_mat,mg_split,restrict,interp] = mg_setup(k,0,op_type,npc,numlev,bc,dim);
    
    %Helmholtz and shifted Laplace matrices
    A    = mg_mat{1};
    Aeps = helmholtz2(k,eps,npf,npf,bc);
    
    %% Sparse FOV of Helmholtz + Shifted Laplace
    tic
    fprintf('beginning computation of lu of Aeps for k=%d  \n', k);
    [L,U] = lu(Aeps);
    UH = L'; LH=U';
    time_lu = toc;
    fprintf('lu factorization for k=%d finished \n', k);
    fprintf('time lu: %f  \n', time_lu);
    
    
    N = length(A);
    Ahat    = @(x) A*(U\(L\x)); %Ahat is the shifted Laplacian
    AhatH   = @(x) UH\(LH\(A'*x));
    H       = @(x) 0.5*(feval(Ahat,x) + feval(AhatH,x)); %Hermitian part of Ahat
    
    
    %% Sparse FOV of deflated shifted Laplacian
    
    %Set up two-level structure
    option.twolevel = true;
    numlev = 2;
    
    
    %Deflation subspace: columns of interpolation operator
    Z = interp{1};
    dim_def = size(Z,2);
    
    %Deflated-shifted Operator PadefA
    Ac = Z'*A*Z; %Coarse operator
    
    %We compute the projection P = I-AZ(Z'AZ)^(-1)Z
    
    fprintf('beginning computation of lu of coarse matrix for k=%d  \n', k);
    tic
    [Lc, Uc] = lu(Ac);
    time_luc = toc;
    
    fprintf('coarse lu factorization for k=%d finished \n', k);
    fprintf('time lu coarse: %f  \n', time_luc);
    
    LcH = Uc'; UcH = Lc';
    
    
   %We compute the projection  P = AZ(Z'AZ)^(-1)AZ
    
    Acinv    =  @(x) Uc\(Lc\x);   %Inverse of coarse Helmholtz
    AcHinv   =  @(x) UcH\(LcH\x); %Inverse of Hermitian transpose of coarse Helmholtz
    
    P  = @(x) A*Z*Acinv(Z'*x);
    PH = @(x) Z*AcHinv(Z'*A'*x);
   
   %Deflated operator
    I        =  speye(N);
    Aepsinv  =  @(x) U\(L\x);     %Inverse of shifted Laplace
    AepsHinv =  @(x) UH\(LH\x);   %Inverse of Hermitian transpose of shifted Laplace
      
    
    %FOV in the Euclidean inner product
    AP  = @(x) Aepsinv(x-feval(P,x));
    APH = @(x) AepsHinv(x) - feval(PH,AepsHinv(x));
    
    fprintf('beginning computation of fov for k=%d \n', k);
    tic
    vmaxH        = rand(N,1);
    [fovAP1,~,~] = parallel_sfov(AP,APH,vmaxH,N,50);
    fovAP1       = 1 + 1i*eps*fovAP1;
    time_fov     = toc;
    fprintf('\n time fov: %f  \n', time_fov);
    
    %close all
    refovAP1 = real(fovAP1); imfovAP1= imag(fovAP1);
    figure(1)
    plot(refovAP1,imfovAP1,'b','LineWidth',3)     %Plot the field of values
    hold on
    plot(0,0,'k+','Markersize',10,'LineWidth',2); 
    axis('equal');
    axis([-0.2 1.2 -0.7 0.7]);
    xlabel('Re(z)','FontSize',14);
    ylabel('Im(z)','FontSize',14);
    set(gca,'Xtick',[ 0   0.5   1],'FontSize',14);
    set(gca,'Ytick',[-0.5 0 0.5 1],'FontSize',14);
    
    
    x = xlabel('$\mathrm{Re}(z)$');                       %x-axis label
    set(x,'Interpreter','latex');
    y = ylabel('$\mathrm{Im}(z)$','interpreter','latex'); %y-axis label
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
