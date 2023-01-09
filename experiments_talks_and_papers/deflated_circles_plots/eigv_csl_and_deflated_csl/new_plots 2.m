%Plots of eigenvalues, ADEF and Deflated CSL
%using new formulas
clear all;
close all;

ppw = 10;  %points per wavelength 

%Parameters for the shifted Laplacian
%We use a shifted Laplacian of the form
% M =  A-i*eps*I, where A is the Helmholtz matrix 
% and eps = factoreps*k^poweps
%
%The standard choice is factoreps = 0.5, poweps=2 
poweps    = 2;
factoreps = 0.5;
%Choose which plots are generated
plot_csl = 'yes';
plot_defcsl = 'yes';
plot_adef = 'no';

% Parameters
kk = [200 500 1000];
%kk = 20;

for i=1:length(kk)
close all
k   = kk(i);

%if pollution = 'no', np = ceil(k^(3/2))
pollution = 'yes';
%otherwise the number of grid points is chosen
%with a fixed number of points per wavelength
np = ceil(ppw*k/(2*pi))-1; 

if strcmp(pollution,'no')
np = ceil(k^(3/2));
end
    
if (mod(np+1,2)==1) 
    np = np+1; 
end
npc = (np-1)/2;

k*(1/np)

%shift for shifted Laplacian
eps = factoreps*k^poweps;

%Eigenvalues of CSL, DCSL,ADEF computed symbolically
[eCSL,eDCSL] = eigSL(k,np,eps);
[eADEFCSL]   = eigADEF1(k,np,eps);

%% Plots: Requires the Matlab2Tikz package!
% stand alone .tex files are saved to the current 
%         MATLAB directory

%Eigenvalues CSL
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
plot(1/2+1/2*cos(t),1/2*sin(t));

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
name1 = strcat('csl_wn',wn,'_ppw',pts, ...
            '_pshift_',powershift,'_fshift_',factorshift,'.tex');
        
        if strcmp(pollution,'no')
            name1 = strcat('csl_wn',wn,'_nopoll', ...
            '_pshift_',powershift,'_fshift_',factorshift,'.tex');
        end

matlab2tikz('filename',name1,'standalone',true); %save .tex file of tikz figure%
end %End of eigenvalues CSL

close all;

if strcmp(plot_defcsl,'yes')

%Eigenvalues of Def CSL
plot(real(eDCSL),imag(eDCSL),'.b','MarkerSize',16); 
axis equal
axis([-1 1.0 -1 1]);
xlabel('Re(z)','FontSize',14);
ylabel('Im(z)','FontSize',14);
set(gca,'Xtick',[-1 -0.5 0 0.5 1],'FontSize',14);
set(gca,'Ytick',[-1 -0.5 0 0.5 1],'FontSize',14);
hold on; 

f = figure;

%setting file names for the .tex files
wn     = num2str(k);  pts = num2str(ppw);
epsshift = num2str(eps);
  
%filename format: wavenumber_pointswavelength_realshift_imagshift.tex
%Tikz Axis formatting
x = xlabel('$\mathrm{Re}(z)$');
% x-axis label
set(x,'Interpreter','latex')
y=ylabel('$\mathrm{Im}(z)$','interpreter','latex'); % x-axis label
set(y,'Interpreter','latex')
figure(1)
name1 = strcat('defcsl','_wn',wn,'_ppw',pts,...
            '_pshift_',powershift,'_fshift_',factorshift,'.tex');

        if strcmp(pollution,'no')
            name1 = strcat('dcsl_wn',wn,'_nopoll', ...
            '_pshift_',powershift,'_fshift_',factorshift,'.tex');
        end
                
matlab2tikz('filename',name1,'standalone',true); %save .tex file of tikz figure%

end %End of plot eigenvalues deflated CSL

close all;

%Begin plot of eigenvalues ADEF+CSL
if strcmp(plot_adef,'yes')

plot(real(eADEFCSL),imag(eADEFCSL),'.b','MarkerSize',16); 
axis equal
axis([-2 2 -2 2]);
xlabel('Re(z)','FontSize',14);
ylabel('Im(z)','FontSize',14);
set(gca,'Xtick',[-2 -1.5 -1 -0.5 0 0.5 1 1.5 2],'FontSize',14);
set(gca,'Ytick',[-2 -1.5 -1 -0.5  0 0.5 1 1.5 2],'FontSize',14);
hold on; 

%setting file names for the .tex files
wn     = num2str(k);  pts = num2str(ppw);
powershift  = num2str(poweps);
factorshift = num2str(factoreps);
  
%filename format: wavenumber_pointswavelength_realshift_imagshift.tex

%Tikz Axis formatting
x = xlabel('$\mathrm{Re}(z)$');
% x-axis label
set(x,'Interpreter','latex')
y=ylabel('$\mathrm{Im}(z)$','interpreter','latex'); % x-axis label
set(y,'Interpreter','latex')
figure(1)
name1 = strcat('adefcsl_wn',wn,'_ppw',pts, ...
            '_pshift_',powershift,'_fshift_',factorshift,'.tex');
        
        if strcmp(pollution,'no')
            name1 = strcat('adefcslcsl_wn',wn,'_nopoll', ...
            '_pshift_',powershift,'_fshift_',factorshift,'.tex');
        end

matlab2tikz('filename',name1,'standalone',true); %save .tex file of tikz figure%
end %End of eigenvalues ADEF + CSL

end
