%Field of values of finite element matrices 1D
%ADEF preconditioner

%% Construction of the matrices
clear all
close all
dim       = 1;
poweps    = 2;
factoreps = 0.5;
bc        = 'som';

%Wavenumber
kmult  =  [10 20 40 80 100];
kk     =  kmult;
iter_SL = zeros(length(kk),2);
minfov  = zeros(length(kk),1);

%Line colors and types for plots
linetyp = {'-','-.',':','--','-'};
%color   = {'r','b','g','k'};
color1  = [0 0 0; 0.5 0 0.5; 0 0 1; 0 0.5 0; 1 0 0];
marker  = {'*','o','.','x'};
opt     = {'m','b','g','k'};
%str2   ={'k','k.-','k--','b*-'};

plot_csl     = 'yes';
plot_gmres   = 'no';
ppw = 0.5;

%if pollution = 'no' the number of points np= ceil(k^(3/2))
pollution = 'no';
fvpts = 60;

%gmres parameters
restart = [];
tol     = 1e-8;
maxit   = 100;

%% Plot of FOV of Shifted Laplace problems
for i=1:length(kk)
    k   = kk(i);
    eps = factoreps*k^poweps;
    
    %choosing the number of points
    if strcmp(pollution,'no')
        npc = ceil(k^2/2);
    end
    
    npf = 2*npc+1;
    
    A    = helmholtzfem(k,npf,0,bc);           %Helmholtz matrix
    Aeps = helmholtzfem(k,npf,eps,bc);         %Shifted Laplace matrix
    M    = mass(npf,bc);
    
    [L,U] = lu(Aeps);
    LH    = U'; UH = L';
    N     = length(A);
    
    %Let C = MAepsinv
    C    = @(x) M*(U\(L\x));
    CH   = @(x) UH\(LH\(M*x));
    HerC = @(x) 0.5*(feval(C,x) + feval(CH,x)); %Hermitian part of C
    
    
    opts.isreal = 0;
    opts.p      = 60;
    opts.issym  = true;
    opts.tol    = 1e-10;
    opts.v0     = rand(N,1);
    fprintf('beginning computation of max eigvalue of H(C) for k = %d  \n', k);
    tic
    [vmaxHerC, eigmaxHerC] =  eigs(HerC,N,1,'LM',opts);
    time_maxeigv = toc;
    fprintf('computation of max eigenvalue of H(C) for k=%d finished \n', k);
    fprintf('time maxeigv: %f  \n', time_maxeigv);
    
    
    %Use the maximum eigenvalue of the Hermitian part of C
    %to initiate computation of fovC
    fprintf('beginning computation of fov for k=%d \n', k);
    tic
    [fovC,~,~] = sfov(C,CH,vmaxHerC,N,60);
    time_fov = toc;
    fprintf('time fov: %f  \n', time_fov);
    
    fovAhat = 1+1i*eps*fovC;
    reFOV  = real(fovAhat); imFOV = imag(fovAhat);
    cvh    = convhull(reFOV,imFOV);
    
    % plotting the fov
    label       = strcat('$k = ', num2str(kmult(i)),'$');
    fovplot(i)  = plot(reFOV,imFOV,'Color',color1(i,:),...
        'LineWidth',2.5,'linestyle',linetyp{i},...
        'DisplayName',label);  
    hold on
    
    %plot(0,0,'k*','Markersize',10,'LineWidth',2);
    %plot(1,0,'b*','Markersize',10,'LineWidth',2);
    axis equal
    axis([-0.2 1.2 -0.7 0.7]);
    xlabel('Re(z)','FontSize',20,'Interpreter','latex');
    ylabel('Im(z)','FontSize',20,'Interpreter','latex');
    h=gca;
 
    
    set(gca,'Xtick',[0 0.5 1],'FontSize',20);
    set(gca,'Ytick',[-0.5 0 0.5],'FontSize',20);
    set(gca,'TickLabelInterpreter', 'tex');
    set(legend,'Interpreter','latex');
    
    
end
L=legend(fovplot);
%set(L,'Interpreter','latex','FontSize',16);

%Filename format:
%1d_fov_csl_kmin_kmax_pointswavelength_realshift_imagshift.tex

kmin  = num2str(min(kmult));
kmax  = num2str(max(kmult));
pts   = num2str(ppw);
powershift  = num2str(poweps);
factorshift = num2str(10*factoreps);

%Tikz Axis formatting
x = xlabel('$\mathrm{Re}(z)$');
%x-axis label
set(x,'Interpreter','latex','fontsize',16)
y=ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',16); % x-axis label
set(y,'Interpreter','latex','fontsize',16)
figure(1)
name1 = strcat('1d_fov_csl_kmin',kmin,'_kmax',kmax,'_ppw',pts, ...
    '_pshift_',powershift,'_fshift_',factorshift,'.tex');

if strcmp(pollution,'no')
    name1 = strcat('1d_fov_csl_kmin',kmin,'_kmax',kmax,'nopoll',pts, ...
        '_pshift_',powershift,'_fshift_',factorshift,'.tex');
end

%save .tex file of tikz figure%
matlab2tikz('filename',name1,'standalone',true,...
    'interpretTickLabelsAsTex',true,...
    'extraaxisoptions',['xlabel style={font=\Large},'...
    'ylabel style={font=\Large},',...
    'legend style={font=\Large},',...
    'ticklabel style={font=\Large}']);

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
