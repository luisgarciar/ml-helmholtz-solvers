%Field of values of finite element matrices 1D
%ADEF preconditioner

%% Construction of the matrices
close all

dim = 1;
poweps    = 2;
factoreps = 10;
bc = 'som';

plot_defcsl = 'no';
%plot_csl_defcsl = 'yes';

%Wavenumber
kmult  =  [10 20 50 80 100];
kk     =  kmult*pi;
pollexp = 1.5;

ppw = 0.5;
%if pollution = 'no' the number of points np= ceil(k^(pollexp))
pollution = 'no';

fvpts = 64;
fovAP = zeros(fvpts, length(kk));


%% Computation of FOV of Deflated Shifted Laplace problems
for i=1:length(kk)
    k   = kk(i);
    eps = factoreps*k^poweps;
    
    %choosing the number of points
    if strcmp(pollution,'no')
        npc = ceil(k^(pollexp)/2);
    end
    
    npf = 2*npc+1;
    
    A    = helmholtzfem(k,npf,0,bc);           %Helmholtz matrix
    Aeps = helmholtzfem(k,npf,eps,bc);         %Shifted Laplace matrix
    
    
    % Sparse FOV of Deflated shifted Laplacian
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
    
    %Deflated operator
    N        = length(A);
    I        = speye(N);
    Aepsinv  =  @(x) U\(L\x);     %Inverse of shifted Laplace
    AepsHinv =  @(x) UH\(LH\x);   %Inverse of Hermitian transpose of shifted Laplace
    Acinv    =  @(x) Uc\(Lc\x);   %Inverse of coarse Helmholtz
    AcHinv   =  @(x) UcH\(LcH\x); %Inverse of Hermitian transpose of coarse Helmholtz
    
    %FOV in the Euclidean inner product
    AP  =  @(x)  M*Aepsinv(x-A*Z*Acinv((Z'*x)));
    APH =  @(x)  AepsHinv(M*x)- Z*AcHinv(Z'*A'*AepsHinv(M*x));
    
    v0 = ones(length(A),1);
    %Field of values of APadef in Euclidean inner product
    [fovAP(:,i),~,~] = sfov(AP,APH,v0,N,64);
    fovAP(:,i) = 1+1i*eps*fovAP(:,i);
    
end


%Line colors and types for plots
linetyp = {'-','-.',':','--','-','-.','-','-.'};
color1  = [0 0 0; 0.5 0 0.5; 0 0 1; 0 0.5 0; 1 0 0; 1 0.5 0; 0 0 0;... 
           0.5 0 0.5];

%% Plotting part
for i=1:length(kk)
    k   = kk(i);   
    eps = factoreps*k^poweps;

    label = strcat('$k = ', num2str(kmult(i)),' \pi$');
    fovplot(i) =  plot(real(fovAP(:,i)),imag(fovAP(:,i)),'Color',color1(i,:),...
                      'LineWidth',2,'linestyle',linetyp{i},...
                      'DisplayName',label);
       
    hold on
    plot(0,0,'k.','Markersize',10,'LineWidth',1);
    plot(1,0,'b.','Markersize',10,'LineWidth',1);
    axis equal
    axis([-0.2 1.2 -0.7 0.7]);
    xlabel('Re(z)','FontSize',16,'Interpreter','latex');
    ylabel('Im(z)','FontSize',16,'Interpreter','latex');
    h=gca;
   
    set(gca,'Xtick',[0 0.5 1],'FontSize',16);
    set(gca,'Ytick',[-0.5 0 0.5],'FontSize',16);
    set(gca,'TickLabelInterpreter', 'tex');

end

L=legend(fovplot);
set(L,'Interpreter','latex','FontSize',16);



%% Saving Figures
%Filename format:
%1d_fov_csl_kmin_kmax_pointswavelength_realshift_imagshift.tex

kmin  = num2str(min(kmult));
kmax  = num2str(max(kmult));
pts   = num2str(ppw);
polexp = num2str(10*pollexp);
powershift  = num2str(poweps);
factorshift = num2str(10*factoreps);

%Tikz Axis formatting
x = xlabel('$\mathrm{Re}(z)$');
%x-axis label
set(x,'Interpreter','latex','fontsize',16)
y=ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',16); % x-axis label
set(y,'Interpreter','latex','fontsize',16)
fig = figure(1);

nametex = strcat('1d_fov_def_csl_kmin',kmin,'_kmax',kmax,'_ppw',pts, ...
    '_pshift_',powershift,'_fshift_',factorshift,'.tex');

if strcmp(pollution,'no')
    nametex = strcat('1d_fov_def_csl_kmin',kmin,'_kmax',kmax,'_nopoll10exp_',polexp, ...
        '_powshift_',powershift,'_10timesfacshift_',factorshift,'.tex');
end

%path for saving the files
currentpath = pwd;
filetex = fullfile(currentpath,'plots','tex',nametex);

%save .tex file of tikz figure%
matlab2tikz('filename',filetex,'standalone',true,...
             'interpretTickLabelsAsTex',true,...
             'extraaxisoptions',['xlabel style={font=\Large},'...
             'ylabel style={font=\Large},',...
             'legend style={font=\Large},',...
             'ticklabel style={font=\Large}']);

nameeps = strcat('1d_fov_def_csl_kmin',kmin,'_kmax',kmax,'_ppw',pts, ...
    '_pshift_',powershift,'_fshift_',factorshift,'.eps');
%nameeps = 'something';

plot_file_eps = fullfile(currentpath,'plots','eps',nameeps);
print(fig,'-depsc',plot_file_eps)         
         
         
         
