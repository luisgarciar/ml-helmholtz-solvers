function plotfovAhat(fovAhat,kmult,opts)
%PLOTFOVAHAT Computes the field of values of 
%preconditioned Helmholtz matrices

%Input:
%kmult: vector of integers 
        %(wavenumbers, to be multiplied by k)
%opts: Data structure with options for the computation
%opts.prec = 'csl' or 'adef': preconditioner
%opts.poweps    = power of the shifted Laplacian
%opts.factoreps = factor of the shifted Laplacian
%opts.dim  = dimension of the problem (1 or 2)
%opts.bc  = Boundary conditions, 'som' or 'dir' 
%opts.disc = Discretization options 
%            'q' for quasioptimal --> n = Ck^2
%            'pf' for pollution free --> n = Ck^1.5
%            'pw' for points per wavelength --> n = Ck  
%
%opts.fvpts = number of points in the field of values
%
%
%Output:
%fovAhat = Matrix of size (opts.fvpts,length(kmult)) containing
%          the boundaries of the field of values of the matrices
%          
%Version 0.1 -- Feb 2019
%
%%

dim = opts.dim;
poweps  = opts.poweps;
factoreps = opts.factoreps;
prec = opts.prec;
disc = opts.disc;

linetyp = {'-','-.',':','--','-','-.','-','-.'};
color1  = [0 0 0; 0.5 0 0.5; 0 0 1; 0 0.5 0; 1 0 0; 1 0.5 0; 0 0 0;... 
           0.5 0 0.5];



%Wavenumber
kk      =  kmult*pi;

switch disc
    case 'q'
        pollexp = 2;
    case 'pf'
        pollexp = 1.5;
    case 'pw'
        pollexp = 1;
end

for i=1:length(kk)
    k   = kk(i);
    eps = factoreps*k^poweps;
   
  
    label = strcat('$k = ', num2str(kmult(i)),' \pi$');
    fovplot(i) =  plot(real(fovAhat(:,i)),imag(fovAhat(:,i)),'Color',color1(i,:),...
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


nametex = strcat('1d_fov',prec,'_kmin',kmin,'_kmax',kmax,'_ppw',pts, ...
    '_pshift_',powershift,'_10timesfacshiftfshift_',factorshift,'.tex');

if strcmp(pollution,'no')
    nametex = strcat('1d_fov',prec,'_kmin',kmin,'_kmax',kmax,'_nopoll10exp_',polexp, ...
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

nameeps = strcat('1d_fov',prec,'_kmin',kmin,'_kmax',kmax,'_ppw',pts, ...
    '_pshift_',powershift,'_10timesfacshift_',factorshift,'.eps');
%nameeps = 'something';

plot_file_eps = fullfile(currentpath,'plots','eps',nameeps);
print(fig,'-depsc',plot_file_eps)      


end