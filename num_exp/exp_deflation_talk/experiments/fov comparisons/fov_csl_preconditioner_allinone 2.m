%Field of values of finite element matrices 1D
%ADEF preconditioner

%% Construction of the matrices
clear all
close all

dim       = 1;
poweps    = 1.5;
factoreps = 1;
bc        = 'som';

%Wavenumber
%kk = 40;
kk      = [20 40 60 80 100];
iter_SL = zeros(length(kk),2);

minfov  = zeros(length(kk),1);
%Line colors and types for plots
linetyp  = {'-','-.',':','--','-'};
linetyp2 = {'-','-','-','-','-'};

color   = {'r','b','g','k'};
color1  = [0 0 0; 0.5 0 0.5; 0 0 1; 0 0.5 0; 1 0 0];
marker  = {'*','o','.','x'};
opt     = {'m','b','g','k'};
%str2={'k','k.-','k--','b*-'};

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
    npf = ceil(ppw*k/(2*pi))-1;
    
    if strcmp(pollution,'no')
        npf = ceil(k^(3/2));
    end
    
    if (mod(npf+1,2)==1)  %set an even number of interior points in 1D
        npf = npf+1;
    end
    npc = (npf-1)/2;
    
    A    = helmholtzfem(k,npf,0,bc);             %Helmholtz matrix
    Aeps = helmholtzfem(k,npf,eps,bc);           %Shifted Laplace matrix
    
   [fovAhat,minfov(i)] = slapfov(A,Aeps,fvpts);   %field of values (complex valued vector)  
   reFOV  = real(fovAhat); imFOV = imag(fovAhat);
   cvh    = convhull(reFOV,imFOV);
   
   %plotting the fov
   label       = strcat('$k = ',num2str(k),'$');
   fovplot(i)  = plot(reFOV(cvh),imFOV(cvh),'Color',color1(i,:),...
                   'LineWidth',2.5,'linestyle',linetyp2{i},...
                  'DisplayName',label);
                
   % fovplot(i)  = plot(reFOV(cvh),imFOV(cvh),'Color',color1(i,:),...
    %                'LineWidth',2,...
     %               'DisplayName',label);        
                
   hold on
   
   plot(0,0,'k*','Markersize',10,'LineWidth',2);
   plot(1,0,'b*','Markersize',10,'LineWidth',2);
   axis equal
   axis([-0.2 1.2 -0.7 0.7]);
   xlabel('Re(z)','FontSize',10,'Interpreter','latex');
   ylabel('Im(z)','FontSize',10,'Interpreter','latex');
   h=gca;
  % [hx,hy] = format_ticks(gca,{'$0$','$0.5$','$1$'},...
   %                      {'$-0.5$','$0$','$0.5$'},...
    %                        [0 0.5 1],[-0.5 0 0.5]);
%       
   
    set(gca,'Xtick',[0 0.5 1],'FontSize',30);
    set(gca,'Ytick',[-0.5 0 0.5],'FontSize',30);
    %set(gca,'TickLabelInterpreter', 'tex');

          
end
L=legend(fovplot);
set(L,'Interpreter','latex','FontSize',40);
          
%Filename format:
%1d_fov_csl_kmin_kmax_pointswavelength_realshift_imagshift.tex
 
kmin  = num2str(min(kk));
kmax  = num2str(max(kk));
pts   = num2str(ppw);
powershift  = num2str(poweps);
factorshift = num2str(10*factoreps);

%Tikz Axis formatting
 x = xlabel('$\mathrm{Re}(z)$');
 %x-axis label
 set(x,'Interpreter','latex','fontsize',30)
 y=ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',30); % x-axis label
 set(y,'Interpreter','latex','fontsize',30)
figure(1)
name1 = strcat('1d_fov_csl_kmin',kmin,'_kmax',kmax,'_ppw',pts, ...
            '_pshift_',powershift,'_fshift_',factorshift,'.tex');
        
if strcmp(pollution,'no')
name1 = strcat('1d_fov_csl_kmin',kmin,'_kmax',kmax,'nopoll', ...
            '_pshift_',powershift,'_fshift_',factorshift,'.tex');
end

%save .tex file of tikz figure%
matlab2tikz('filename',name1,'standalone',true,...
   'interpretTickLabelsAsTex',true,...
      'extraaxisoptions',['xlabel style={font=\LARGE},'...
                       'ylabel style={font=\LARGE},',...
                       'legend style={font=\LARGE},',...
                       'ticklabel style={font=\HUGE}']);
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
   


