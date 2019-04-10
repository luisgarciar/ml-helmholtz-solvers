function plotfovAhat(fovAhat,kmult,opts)
%PLOTFOVAHAT Computes the field of values of 
%preconditioned Helmholtz matrices

%Input:
%kmult: vector of integers 
        %(wavenumbers, to be multiplied by pi)
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
    xlabel('Re(z)','FontSize',14,'Interpreter','latex');
    ylabel('Im(z)','FontSize',14,'Interpreter','latex');
    h=gca;
   
    set(gca,'Xtick',[0 0.5 1],'FontSize',14);
    set(gca,'Ytick',[-0.5 0 0.5],'FontSize',14);
    set(gca,'TickLabelInterpreter', 'tex');
 
end

L=legend(fovplot);
set(L,'Interpreter','latex','FontSize',14);

kmin        = num2str(min(kmult));
kmax        = num2str(max(kmult));
polexp      = num2str(10*pollexp);
powershift  = num2str(poweps);
factorshift = num2str(10*factoreps);

%Tikz Axis formatting
x = xlabel('$\mathrm{Re}(z)$');
%x-axis label
set(x,'Interpreter','latex','fontsize',14)
y=ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',14); % x-axis label
set(y,'Interpreter','latex','fontsize',14)
fig = figure(1);

nametex = strcat('1dfv',prec,'kmin',kmin,'kmax',kmax,'d',disc,'px',polexp, ...
       'pws',powershift,'fs',factorshift,'.tex');

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

nameeps = strcat('1dfv',prec,'kmin',kmin,'kmax',kmax,'d',disc,'px',polexp, ...
       'pws',powershift,'fs',factorshift);
   
plot_file_eps = fullfile(currentpath,'plots','eps',nameeps);
print(fig,'-depsc',plot_file_eps)      


end