%Field of values of finite element matrices 1D
%ADEF preconditioner

%% Construction of the matrices
clear all
close all

dim = 1;
poweps    = 2;
factoreps = 10;
bc = 'som';

plot_defcsl = 'no';
plot_csl_defcsl = 'yes';

%Wavenumber
%kk      = [20 40 60 80 100];
kk = 200;

%Colors for plots
minfov  = zeros(length(kk),1);
linetyp = {'-','-','-.',':','-'};
%color   = {'r','b','g','k'};
color1  = [0 0 0; 0.5 0 0.5; 0 0 1; 0 0.5 0; 1 0 0];
marker  = {'*','o','.','x'};
opt     = {'m','b','g','k'};

ppw = 0.5;
%if pollution = 'no' the number of points np= ceil(k^(3/2))
pollution = 'no';

fvpts = 64;


%% Plot of FOV of Shifted Laplace problems
for i=1:length(kk)
    k   = kk(i);
    eps = factoreps*k^poweps;
    
    %choosing the number of points
    npf = ceil(ppw*k/(2*pi))-1;
    if strcmp(pollution,'no')
        npf = ceil(k^(2));
    end
    
    if (mod(npf+1,2)==1)  %set an even number of interior points in 1D
        npf = npf+1;
    end
    
    npc  = (npf-1)/2;
    A    = helmholtzfem(k,npf,0,bc);           %Helmholtz matrix
    Aeps = helmholtzfem(k,npf,eps,bc);         %Shifted Laplace matrix
    
    
    %dim  = 2;
    %A    = helmholtz2(k,0,npf,npf,bc);
    %Aeps = helmholtz2(k,eps,npf,npf,bc);
    
    %% Sparse FOV of Deflated shifted Laplacian
    R  = fwrestrictionfem(npf,dim,bc);
    Z  = R';             %Prolongation. Deflation subspace: columns of Z
    dim_def = size(Z,2);
    
    [L,U] = lu(Aeps);
    LH = U'; UH = L';
    
    %Deflated-shifted Operator PadefA
    M = mass(npf,bc);
    Ac = Z'*A*Z; %Coarse operator
    [Lc, Uc] = lu(Ac);
    LcH = Uc'; UcH = Lc';
    sqrtM = sqrtm(full(M));
    
    
    %Deflated operator
    N        = length(A);
    I        = speye(N);
    Aepsinv  =  @(x) U\(L\x);     %Inverse of shifted Laplace
    AepsHinv =  @(x) UH\(LH\x);   %Inverse of Hermitian transpose of shifted Laplace
    Acinv    =  @(x) Uc\(Lc\x);   %Inverse of coarse Helmholtz
    AcHinv   =  @(x) UcH\(LcH\x); %Inverse of Hermitian transpose of coarse Helmholtz
    
    %FOV in the Minv inner product:
    %AP  =  @(x)  sqrtM*Aepsinv((sqrtM*x-A*Z*Acinv((Z'*sqrtM*x))));
    %APH =  @(x)  sqrtM*(AepsHinv(sqrtM*x)- Z*AcHinv(Z'*A'*AepsHinv(sqrtM*x)));
    
    %FOV in the Euclidean inner product
    AP  =  @(x)  M*Aepsinv(x-A*Z*Acinv((Z'*x)));
    APH =  @(x)  AepsHinv(M*x)- Z*AcHinv(Z'*A'*AepsHinv(M*x));
    
    v0 = ones(length(A),1);
    %Field of values of APadef in Euclidean inner product
    [fovAP,~,~] = sfov(AP,APH,v0,N,64);
    fovAP = 1+1i*eps*fovAP;
    
    close all
    if strcmp(plot_defcsl,'yes')
        label   = strcat('$k = ',num2str(k),'$');
        plot(real(fovAP),imag(fovAP),'Color',color1(i,:),...
            'LineWidth',4,'linestyle',linetyp{i},...
            'DisplayName',label);
        hold on
        plot(0,0,'k*','Markersize',10,'LineWidth',2);
        plot(1,0,'b*','Markersize',10,'LineWidth',2);
        axis equal
        axis([-0.2 1.2 -0.7 0.7]);

       % axis([0.7 1.3 -0.3 0.3]);
        xlabel('Re(z)','FontSize',30,'Interpreter','latex');
        ylabel('Im(z)','FontSize',30,'Interpreter','latex');
        h=gca;
        
        set(gca,'Xtick',[0 0.5 1],'FontSize',30);
        set(gca,'Ytick',[-0.5 0 0.5],'FontSize',30);
        %set(gca,'Xtick',[0.7 1 1.3],'FontSize',30);
        %set(gca,'Ytick',[-0.3 0 0.4],'FontSize',30);
        
        kmin  = num2str(min(kk));
        kmax  = num2str(max(kk));
        pts   = num2str(ppw);
        powershift  = num2str(poweps);
        factorshift = num2str(10*factoreps);
         
        wn          = num2str(k);  pts = num2str(ppw);
        powershift  = num2str(poweps);
        factorshift = num2str(10*factoreps);
        
        x = xlabel('$\mathrm{Re}(z)$');
        %         % x-axis label
        set(x,'Interpreter','latex')
        y=ylabel('$\mathrm{Im}(z)$','interpreter','latex'); % x-axis label
        set(y,'Interpreter','latex')
        figure(1)
        name1 = strcat('1d_fov_dcsl_kmin',kmin,'_kmax',kmax,'_ppw',pts, ...
            '_pshift_',powershift,'_fshift_',factorshift,'.tex');
        
        if strcmp(pollution,'no')
            name1 = strcat('1d_fov_dcsl_wn',wn,'_nopoll', ...
                '_pshift_',powershift,'_fshift_',factorshift,'.tex');
        end
        %
        matlab2tikz('filename',name1,'standalone',true,...
            'interpretTickLabelsAsTex',true,...
            'extraaxisoptions',['xlabel style={font=\LARGE},'...
            'ylabel style={font=\LARGE},',...
            'legend style={font=\LARGE},',...
            'ticklabel style={font=\HUGE}']);    
    end %End of FOV DCSL
    
    if strcmp(plot_csl_defcsl,'yes')
        
        [fovAhat,minfov(i)] = slapfov(A,Aeps,fvpts);   %field of values (complex valued vector)
        reFOV  = real(fovAhat); imFOV = imag(fovAhat);
        cvh    = convhull(reFOV,imFOV);
        
        close all
        %plotting the fov
        label       = strcat('$k = ',num2str(k),'$');
        plot(reFOV(cvh),imFOV(cvh),'Color',color1(i,:),...
            'LineWidth',2.5,'linestyle',linetyp{i},...
            'DisplayName',label);
        
        hold on
        plot(real(fovAP),imag(fovAP),'Color',color1(i+1,:),...
            'LineWidth',2.5,'linestyle',linetyp{i+1},...
            'DisplayName',label);
        plot(0,0,'k*','Markersize',10,'LineWidth',2);
        plot(1,0,'b*','Markersize',10,'LineWidth',2);
        axis equal
        axis([-0.2 1.2 -0.7 0.7]);
        xlabel('Re(z)','FontSize',30,'Interpreter','latex');
        ylabel('Im(z)','FontSize',30,'Interpreter','latex');
        h=gca;
        
        set(gca,'Xtick',[0 0.5 1],'FontSize',30);
        set(gca,'Ytick',[-0.5 0 0.5],'FontSize',30);
        %set(gca,'TickLabelInterpreter', 'tex');
        
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
        
    matlab2tikz('filename',name1,'standalone',true,...
            'interpretTickLabelsAsTex',true,...
            'extraaxisoptions',['xlabel style={font=\LARGE},'...
            'ylabel style={font=\LARGE},',...
            'legend style={font=\LARGE},',...
            'ticklabel style={font=\HUGE}']);    
    end %End of FOV CSL & DCSL
    %close all;
end
