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
%kk = 20;
kk = [140 180];


plot_csl    = 'yes';
plot_defcsl = 'yes';
plot_csl_defcsl = 'yes';
ppw = 0.5;
%if pollution = 'no' the number of points np= ceil(k^(3/2))
pollution = 'no';

for i=1:length(kk)
    close all
    
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
    
    A  = helmholtzfem(k,npf,0,bc); %Helmholtz matrix
    Aeps = helmholtzfem(k,npf,eps,bc); %Shifted Laplace matrix
    
    %% Sparse FOV of Helmholtz + Shifted Laplace
    [L,U] = lu(Aeps);
    LH = U'; UH=L';
    N=length(A);
    
    Ahat     = @(x) A*(U\(L\x));
    AhatH    = @(x) UH\(LH\(A'*x));
    H        = @(x) 0.5*(feval(Ahat,x) + feval(AhatH,x)); %Hermitian part of Ahat
    
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
    
    if strcmp(plot_csl,'yes')
        plot(reFOV(cvh), imFOV(cvh),'b','LineWidth',4)      % Plot the field of values
        hold on
        plot(0,0,'k+','Markersize',10,'LineWidth',4)
        axis equal
        axis([-0.2 1.2 -0.7 0.7]);
        xlabel('Re(z)','FontSize',14);
        ylabel('Im(z)','FontSize',14);
        set(gca,'Xtick',[0 0.5 1],'FontSize',14);
        set(gca,'Ytick',[-0.5 0 0.5 1],'FontSize',14);
        
        wn     = num2str(k);  pts = num2str(ppw);
        powershift  = num2str(poweps);
        factorshift = num2str(10*factoreps);
        
        %Filename format:
        %wavenumber_pointswavelength_realshift_imagshift.tex
        
        %Tikz Axis formatting
        x = xlabel('$\mathrm{Re}(z)$');
        % x-axis label
        set(x,'Interpreter','latex')
        y=ylabel('$\mathrm{Im}(z)$','interpreter','latex'); % x-axis label
        set(y,'Interpreter','latex')
        figure(1)
        name1 = strcat('1d_fov_csl_wn',wn,'_ppw',pts, ...
            '_pshift_',powershift,'_fshift_',factorshift,'.tex');
        
        if strcmp(pollution,'no')
            name1 = strcat('1d_fov_csl_wn',wn,'_nopoll', ...
                '_pshift_',powershift,'_fshift_',factorshift,'.tex');
        end
        
        matlab2tikz('filename',name1,'standalone',true); %save .tex file of tikz figure%
    end %End of FOV CSL
    close all;
    
    
    %% Sparse FOV of Deflated shifted Laplacian
    R  = fwrestrictionfem(npf,dim,bc);
    Z  = R';             %Prolongation. Deflation subspace: columns of Z
    dim_def = size(Z,2);
    
    %Deflated-shifted Operator PadefA
    M = mass(npf,bc);
    Ac = Z'*A*Z; %Coarse operator
    [Lc, Uc] = lu(Ac);
    LcH = Uc'; UcH = Lc';
    sqrtM = sqrtm(full(M));
    
    %Deflated operator
    I     = speye(N);
    Aepsinv  =  @(x) U\(L\x);     %Inverse of shifted Laplace
    AepsHinv =  @(x) UH\(LH\x);   %Inverse of Hermitian transpose of shifted Laplace
    Acinv    =  @(x) Uc\(Lc\x);   %Inverse of coarse Helmholtz
    AcHinv   =  @(x) UcH\(LcH\x); %Inverse of Hermitian transpose of coarse Helmholtz
    
    AP  =    @(x)  sqrtM*Aepsinv((sqrtM*x-A*Z*Acinv((Z'*sqrtM*x))));
    APH =    @(x)  sqrtM*(AepsHinv(sqrtM*x)- Z*AcHinv(Z'*A'*AepsHinv(sqrtM*x)));
    
    %Field of values of APadef in Minv inner product
    [fovAP,~,~] = sfov(AP,APH,vmaxH,N,50);
    fovAP = 1+1i*eps*fovAP;
    
    close all
    if strcmp(plot_defcsl,'yes')
        close all
        refovAP = real(fovAP); imfovAP= imag(fovAP);
        plot(refovAP,imfovAP,'k','LineWidth',4)      % Plot the field of values
        hold on 
        hold on
        plot(0,0,'k+','Markersize',10,'LineWidth',4)
        axis('equal');
        axis([-0.2 1.2 -0.7 0.7]);
        xlabel('Re(z)','FontSize',14);
        ylabel('Im(z)','FontSize',14);
        set(gca,'Xtick',[0 0.5 1],'FontSize',14);
        set(gca,'Ytick',[-0.5 0 0.5 1],'FontSize',14);
        
        x = xlabel('$\mathrm{Re}(z)$');
        % x-axis label
        set(x,'Interpreter','latex')
        y=ylabel('$\mathrm{Im}(z)$','interpreter','latex'); % x-axis label
        set(y,'Interpreter','latex')
        figure(1)
        name1 = strcat('1d_fov_dcsl_wn',wn,'_ppw',pts, ...
            '_pshift_',powershift,'_fshift_',factorshift,'.tex');
        
        if strcmp(pollution,'no')
            name1 = strcat('1d_fov_dcsl_wn',wn,'_nopoll', ...
                '_pshift_',powershift,'_fshift_',factorshift,'.tex');
        end
        
        matlab2tikz('filename',name1,'standalone',true); %save .tex file of tikz figure%
    end %End of FOV DCSL
    
    if strcmp(plot_csl_defcsl,'yes')
        close all
        plot(reFOV(cvh), imFOV(cvh),'b','LineWidth',4);      % Plot the field of values
        hold on
        plot(0,0,'k+','Markersize',10,'LineWidth',4)
        hold on
        plot(refovAP,imfovAP,'k','LineWidth',4);      % Plot the field of values
        axis equal
        axis([-0.2 1.2 -0.7 0.7]);
        xlabel('Re(z)','FontSize',14);
        ylabel('Im(z)','FontSize',14);
        set(gca,'Xtick',[0 0.5 1],'FontSize',14);
        set(gca,'Ytick',[-0.5 0 0.5 1],'FontSize',14);
        
        wn     = num2str(k);  pts = num2str(ppw);
        powershift  = num2str(poweps);
        factorshift = num2str(10*factoreps);
        
        %Filename format:
        %wavenumber_pointswavelength_realshift_imagshift.tex
        
        %Tikz Axis formatting
        x = xlabel('$\mathrm{Re}(z)$');
        % x-axis label
        set(x,'Interpreter','latex')
        y=ylabel('$\mathrm{Im}(z)$','interpreter','latex'); % x-axis label
        set(y,'Interpreter','latex')
        figure(1)
        name1 = strcat('1d_fov_csl_dcsl_wn',wn,'_ppw',pts, ...
            '_pshift_',powershift,'_fshift_',factorshift,'.tex');
        
        if strcmp(pollution,'no')
            name1 = strcat('1d_fov_csl_dcsl_wn',wn,'_nopoll', ...
                '_pshift_',powershift,'_fshift_',factorshift,'.tex');
        end
        
        matlab2tikz('filename',name1,'standalone',true); %save .tex file of tikz figure%
    end %End of FOV CSL & DCSL
    close all;
 
    
end







%% FOV of Helmholtz + Shifted Laplace with full matrices
% Ahat = full(Aeps\A);
% N    = length(A);
%
% %Field of values of Ahat Sl-preconditioned Helmholtz matrix
% [fvAhat, eigAhat] = fv(Ahat,1,32,1);
%
% figure(1)
% plot(real(fvAhat), imag(fvAhat),'b') %Plot the FOV of Padef
% hold on
% plot(real(eigAhat), imag(eigAhat),'r+')    %Plot the eigenvalues too.
% %axis('equal');

%Old Stuff
%% FOV of Deflated Shifted Laplace with full matrices
%
% R = fwrestrictionfem(npf,dim,bc);
% Z = R';             %Prolongation. Deflation subspace: columns of Z
% dim_def = size(Z,2);
%
% %Deflated-shifted Operator PadefA
% M = mass(npf,bc);
% sqrtM = sqrtm(full(M));
% I = eye(length(A));
% P = full(Aeps\(I-A*Z*((Z'*A*Z)\Z')));
% P = sqrtM*P*sqrtM;
%
% %Field of values of PadefA in Minv inner product
% [fvAP] = 1+1i*eps*fv(P,1,32,1);
%
% %Plots
% figure(2)
% plot(real(fvAhat), imag(fvAhat),'b') %Plot the FOV of Padef
% hold on
% plot(real(eigAhat), imag(eigAhat),'r+')    %Plot the eigenvalues too.
% axis('equal');
% plot(0,0,'+k')
% plot(real(fvAP), imag(fvAP),'k*') %Plot the FOV of Padef





%%
%hold on
%plot(real(eigAP), imag(eigAP),'r+')    %Plot the eigenvalues too.

%
% %% 1D Dirichlet Example
% dim = 1;
% npc = 4;
% bc = 'dir';
% %wavenumber and imaginary shift of shifted Laplacian
% k   = 200;  eps = 0.5*k^2; %Helmholtz problem
% ppw = 15;   %number of points per wavelength
% [npf,numlev] = fd_npc_to_npf(npc,k,ppw);  %number of points in finest grid (1D)
%
% A = helmholtz(k,0,npf,bc);
% M = helmholtz(k,eps,npf,bc);  %Shifted Laplace matrix
%
% R = fwrestriction(npf,dim,bc);
% Z = R'; %Prolongation. Deflation subspace: columns of Z
% dim_def = size(Z,2);
%
% % Field of values of restricted matrix
% %To compute the field of values of the restricted matrix
% %(ADEF) we need a basis Y for the orthogonal complement of
% %the columns of Z, i.e., the nullspace of R.
%
% Y = null(full(R));
%
% %Restricted matrix
% Ainvc = (Y'*(A\Y));
% %Restricted ADEF matrix (only upper block)
% B     = (Y'*(M\Y))*((Y'*(A\Y))\speye(length(Ainvc)));
% [FV_B, eigB] = fv(B,1,32);
%
% %Explicit Deflation operator
% Q = Z*((Z'*A*Z)\Z');
% %Deflated-shifted Operator PadefA
% PadefA =  M\(A-A*Z*((Z'*A*Z)\Z'*A))+ Z*((Z'*A*Z)\Z'*A);
%
% [U,R] = qr(Z,0);
% Porth = U*U';
%
% I = eye(length(A));
% %PadefAcorr = PadefA-Porth*PadefA*(I-Porth);
% PadefAcorr = (I-Porth)*PadefA + Porth;
%
% %Plots
% %Field of values of PadefA
% [FV_PadefA, eigPadefA] = fv(full(PadefA),1,32,1);
% [FV_PadefAcorr, eigPadefAcorr] = fv(full(PadefAcorr),1,32,1);
%
%
% plot(real(FV_PadefA), imag(FV_PadefA),'b') %Plot the FOV of Padef
% %plot(real(eigPadefA), imag(eigPadefA), 'r+')    %Plot the eigenvalues too.
% axis('equal');
% hold on
% %plot(real(FV_PadefAcorr), imag(FV_PadefAcorr),'k')   %Plot the FOV of Padefcorr
% %plot(real(eigPadefAcorr), imag(eigPadefAcorr), 'kx') %Plot the eigenvalues too.
% plot(real(FV_B),imag(FV_B),'r')     %Plot the field of values of B
% %plot(real(eigB),imag(eigB), 'k.')  %Plot the eigenvalues too.
% %plot(0,0,'+k')

%
% %% 2D Sommerfeld Example
% clear all
% dim = 2;
% npc = 4;
% bc = 'som';
% %wavenumber and imaginary shift of shifted Laplacian
% k   = 10;  eps = 0.5*k^2; %Helmholtz problem
% ppw = 10;   %number of points per wavelength
% [npf,numlev] = fd_npc_to_npf(npc,k,ppw);  %number of points in finest grid (1D)
%
% A = helmholtz2(k,0,npf,npf,bc);
% M = helmholtz2(k,eps,npf,npf,bc);  %Shifted Laplace matrix
%
% R = fwrestriction(npf,dim,bc);
% Z = R'; %Prolongation. Deflation subspace: columns of Z
% dim_def = size(Z,2);
%
% % Field of values of restricted matrix
% %To compute the field of values of the restricted matrix
% %(ADEF) we need a basis Y for the orthogonal complement of
% %the columns of Z, i.e., the nullspace of R.
%
% Y = null(full(R));
%
% %Restricted matrix
% Ainvc = (Y'*(A\Y));
% %Restricted ADEF matrix (only upper block)
% B     = (Y'*(M\Y))*((Y'*(A\Y))\speye(length(Ainvc)));
% [FV_B, eigB] = fv(B,1,32);
%
% %Explicit Deflation operator
% Q = Z*((Z'*A*Z)\Z');
% %Deflated-shifted Operator PadefA
% PadefA =  M\(A-A*Z*((Z'*A*Z)\Z'*A))+ Z*((Z'*A*Z)\Z'*A);
%
% [U,R] = qr(Z,0);
% Porth = U*U';
%
% I = eye(length(A));
% %PadefAcorr = PadefA-Porth*PadefA*(I-Porth);
% PadefAcorr = (I-Porth)*PadefA + Porth;
%
% %Plots
% %Field of values of PadefA
% [FV_PadefA, eigPadefA] = fv(full(PadefA),1,32,1);
% [FV_PadefAcorr, eigPadefAcorr] = fv(full(PadefAcorr),1,32,1);
%
%
% plot(real(FV_PadefA), imag(FV_PadefA),'b') %Plot the FOV of Padef
% hold on
% plot(real(eigPadefA), imag(eigPadefA), 'r+')    %Plot the eigenvalues too.
% axis('equal');
% plot(real(FV_PadefAcorr), imag(FV_PadefAcorr),'k')   %Plot the FOV of Padefcorr
% plot(real(eigPadefAcorr), imag(eigPadefAcorr), 'kx') %Plot the eigenvalues too.
% plot(real(FV_B),imag(FV_B),'gx')     %Plot the field of values of B
% plot(real(eigB),imag(eigB), 'k.')  %Plot the eigenvalues too.
% plot(0,0,'+k')
%
%
% %% Testing if there is a difference between the matrices in GMRES
% tol   = 1e-10;
% maxit = length(A);
% b     = randn(length(A),1);
%
% [x1,flag1,relres1,iter1,resvec1] = gmres(PadefA,b,[],tol,maxit,[]);
% [x2,flag2,relres2,iter2,resvec2] = gmres(PadefAcorr,b,[],tol,maxit,[]);
%
% iter1
% iter2
%
% figure(3)
% semilogy(1:(iter2(2)+1),resvec2'/resvec2(1),'r-+')
% hold on
% semilogy(1:(iter1(2)+1),resvec1'/resvec1(1),'b-+')