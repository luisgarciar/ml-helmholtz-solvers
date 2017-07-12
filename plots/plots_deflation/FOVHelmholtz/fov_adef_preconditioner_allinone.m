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
kk = [15 20 30 60];
linetyp={'-','.','--','-'};
color = {'m','b','g','k'};
marker = {'*','o','.','x'};
opt = {'m','b','g','k'};
%str2={'k','k.-','k--','b*-'};

plot_csl    = 'yes';
plot_defcsl = 'yes';
plot_csl_defcsl = 'yes';
ppw = 0.5;
%if pollution = 'no' the number of points np= ceil(k^(3/2))
pollution = 'no';

fvpts = 50;


%% Plot of FOV of Shifted Laplace problems
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
    
    A  = helmholtzfem(k,npf,0,bc); %Helmholtz matrix
    Aeps = helmholtzfem(k,npf,eps,bc); %Shifted Laplace matrix
    
    % Sparse FOV of Helmholtz + Shifted Laplace
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
    [fovAhat,~,~] = sfov(Ahat,AhatH,vmaxH,N,fvpts);
    
    %Find the convex hull of the FOV and plot it
     reFOV = real(fovAhat); imFOV = imag(fovAhat);
        cvh = convhull(reFOV,imFOV);
    
    if strcmp(plot_csl,'yes')
        lt = char(linetyp(i));
        cl = char(color(i));
        mk = char(marker(i));
        op = char(opt(i));
        
        %plot(reFOV(cvh),imFOV(cvh),'Color',cl,'Linestyle',lt,'LineWidth',2,'Markersize',2)      % Plot the field of values
        plot(reFOV(cvh),imFOV(cvh),op,'LineWidth',2)
        hold on
        plot(0,0,'k+','Markersize',6,'LineWidth',2)
        axis equal
        axis([-0.2 1.2 -0.7 0.7]);
        xlabel('Re(z)','FontSize',14);
        ylabel('Im(z)','FontSize',14);
        set(gca,'Xtick',[0 0.5 1],'FontSize',14);
        set(gca,'Ytick',[-0.5 0 0.5 1],'FontSize',14);
        
%         wn     = num2str(k);  pts = num2str(ppw);
%         powershift  = num2str(poweps);
%         factorshift = num2str(10*factoreps);
%         
%         %Filename format:
%         %wavenumber_pointswavelength_realshift_imagshift.tex
%         
%         %Tikz Axis formatting
%         x = xlabel('$\mathrm{Re}(z)$');
%         % x-axis label
%         set(x,'Interpreter','latex')
%         y=ylabel('$\mathrm{Im}(z)$','interpreter','latex'); % x-axis label
%         set(y,'Interpreter','latex')
%         figure(1)
%         name1 = strcat('1d_fov_csl_wn',wn,'_ppw',pts, ...
%             '_pshift_',powershift,'_fshift_',factorshift,'.tex');
%         
%         if strcmp(pollution,'no')
%             name1 = strcat('1d_fov_csl_wn',wn,'_nopoll', ...
%                 '_pshift_',powershift,'_fshift_',factorshift,'.tex');
%         end
%         
%         matlab2tikz('filename',name1,'standalone',true); %save .tex file of tikz figure%
    end %End of FOV CSL
    
end

%%

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
    [fovAP,~,~] = sfov(AP,APH,vmaxH,N,64);
    fovAP = 1+1i*eps*fovAP;
    
    close all
    if strcmp(plot_defcsl,'yes')
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
        
%         x = xlabel('$\mathrm{Re}(z)$');
%         % x-axis label
%         set(x,'Interpreter','latex')
%         y=ylabel('$\mathrm{Im}(z)$','interpreter','latex'); % x-axis label
%         set(y,'Interpreter','latex')
%         figure(1)
%         name1 = strcat('1d_fov_dcsl_wn',wn,'_ppw',pts, ...
%             '_pshift_',powershift,'_fshift_',factorshift,'.tex');
%         
%         if strcmp(pollution,'no')
%             name1 = strcat('1d_fov_dcsl_wn',wn,'_nopoll', ...
%                 '_pshift_',powershift,'_fshift_',factorshift,'.tex');
%         end
%         
%         matlab2tikz('filename',name1,'standalone',true); %save .tex file of tikz figure%
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
 
    



