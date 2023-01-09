clear all;
close all;

ppw = 10;  %points per wavelength 

%Parameters for the shifted Laplacian
%We use a shifted Laplacian of the form
% M =  A-i*eps*I, where A is the Helmholtz matrix 
% and eps = factoreps*k^poweps
%The standard choice is factoreps = 0.5, poweps=2 
poweps    = 2;
factoreps = 0.5;


%Choose which plots are generated
plot_csl    = 'yes';
plot_defcsl = 'yes';

% Parameters


kk = [200 500 1000];

for i=1:length(kk)
close all
k   = kk(i);
eps = factoreps*(k^poweps);

%if pollution = 'no', np = ceil(k^(3/2))
pollution = 'yes';
%otherwise the number of grid points is chosen
%with a fixed number of points per wavelength
np = ceil(ppw*k/(2*pi))-1; 

if strcmp(pollution,'no')
np = ceil(k^(3/2));
end
    
if (mod(np,2)==1) 
    np = np+1; 
end
npc = (np-1)/2;

dim = 1;

A = helmholtzfem(k,np,0); %Helmholtz matrix
M = helmholtzfem(k,np,eps); %Shifted Laplace matrix

Ahat = M\A;
eCSL = eig(full(Ahat));

if strcmp(plot_csl,'yes')
    plot(real(eCSL),imag(eCSL),'.b','MarkerSize',16); 
    axis equal
    axis([-1 1.0 -1 1]);
    xlabel('Re(z)','FontSize',14);
    ylabel('Im(z)','FontSize',14);
    set(gca,'Xtick',[-1 -0.5 0 0.5 1],'FontSize',14);
    set(gca,'Ytick',[-1 -0.5 0 0.5 1],'FontSize',14);
    hold on; 
    t=linspace(0,2*pi,100);
    plot(1/2+1/2*cos(t),1/2*sin(t),'k');

%setting file names for the .tex files
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
name1 = strcat('cslfem_wn',wn,'_ppw',pts, ...
            '_pshift_',powershift,'_fshift_',factorshift,'.tex');
        
        if strcmp(pollution,'no')
            name1 = strcat('cslfem_wn',wn,'_nopoll', ...
            '_pshift_',powershift,'_fshift_',factorshift,'.tex');
        end

matlab2tikz('filename',name1,'standalone',true); %save .tex file of tikz figure%
end %End of eigenvalues dCSL

close all;

R = fwrestrictionfem(np,dim);
Z = R'; %Prolongation. Deflation subspace: columns of Z
dim_def = size(Z,2);

%Field of values of restricted matrix 
%To compute the field of values of the restricted matrix
%(DEFCSL) we need a basis Y for the orthogonal complement of
%the columns of Z, i.e., the nullspace of R.

Y = null(full(R));
%Restricted matrix
B = (Y'*(Ahat\Y));
[FV_B, eigB] = fv(B,1,32,1);

if strcmp(plot_defcsl,'yes')
    plot(real(eigB),imag(eigB),'.b','MarkerSize',16); 
    axis equal
    axis([-1 1.0 -1 1]);
    xlabel('Re(z)','FontSize',14);
    ylabel('Im(z)','FontSize',14);
    set(gca,'Xtick',[-1 -0.5 0 0.5 1],'FontSize',14);
    set(gca,'Ytick',[-1 -0.5 0 0.5 1],'FontSize',14);
    hold on; 
    t=linspace(0,2*pi,100);
    plot(1/2+1/2*cos(t),1/2*sin(t),'k');
       
    %setting file names for the .tex files
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
    name1 = strcat('dcslfem_wn',wn,'_ppw',pts, ...
                '_pshift_',powershift,'_fshift_',factorshift,'.tex');

        if strcmp(pollution,'no')
            name1 = strcat('dcslfem_wn',wn,'_nopoll', ...
            '_pshift_',powershift,'_fshift_',factorshift,'.tex');
        end

matlab2tikz('filename',name1,'standalone',true); %save .tex file of tikz figure%
    
end    
    
end


