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
kk = 40;
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
    
    if (mod(npf+1,2)==1)  %set an even number of interior points in 1D
        npf = npf+1;
    end
    npc = (npf-1)/2;
    
    A    = helmholtz2(k,0,npf,npf,bc);           %Helmholtz matrix
    Aeps = helmholtz2(k,eps,npf,npf,bc);         %Shifted Laplace matrix
    
   [fovAhat,minfov(i)] = slapfov(A,Aeps,fvpts);   %field of values (complex valued vector)  
   reFOV  = real(fovAhat); imFOV = imag(fovAhat);
   cvh    = convhull(reFOV,imFOV);
   
   
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
   


