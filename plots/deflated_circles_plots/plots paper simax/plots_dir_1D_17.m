% Plots for SIMAX paper 2017
clear all;
close all;

% Parameters
kk = [100 200 500];
mineigvCSL = zeros(length(kk),1);
mineigvDCSL = zeros(length(kk),1);

for i=1:length(kk)
close all
k   = kk(i);

pollution = 'yes'; %if pollution = 'no', n = k^(3/2)

ppw       = 10;
poweps    = 2;
factoreps = 0.5;
eps = factoreps*k^poweps;
np  = ceil(ppw*k/(2*pi))-1;
dim = 1;
bc  = 'dir';

if strcmp(pollution,'no')
np = ceil(k^(3/2));
end
    
if (mod(np+1,2)==1) 
    np = np+1; 
end
npc = (np-1)/2;

%Eigenvalues of CSL and DCSL computed symbolically
[eCSL,eDCSL] = eigSL(k,np,eps); 

mineigvCSL(i) =  abs(min(eCSL));
mineigvDCSL(i) = abs(min(eDCSL));

%% Plots for paper
%f = figure; 

%Eigenvalues CSL
plot(real(eCSL),imag(eCSL),'.b','MarkerSize',16); 
axis equal
axis([-1 1.0 -1 1]);
xlabel('Re(z)','FontSize',14);
ylabel('Im(z)','FontSize',14);
set(gca,'Xtick',[-1 -0.5 0 0.5 1],'FontSize',14);
set(gca,'Ytick',[-1 -0.5 0 0.5 1],'FontSize',14);
hold on; 
hold on; 
t=linspace(0,2*pi,100);
plot(1/2+1/2*cos(t),1/2*sin(t),'k');

f = figure; 
%setting file names for the eps images
wn     = num2str(k);  pts = num2str(ppw);
powershift  = num2str(poweps);
factorshift = num2str(10*factoreps);
%epsshift = num2str(eps);
  
%filename format: wavenumber_pointswavelength_realshift_imagshift.tex
%plot(real(eCSL),imag(eCSL),'b.')

%Tikz Axis formatting
x = xlabel('$\mathrm{Re}(z)$');
% x-axis label
set(x,'Interpreter','latex')
y=ylabel('$\mathrm{Im}(z)$','interpreter','latex'); % x-axis label
set(y,'Interpreter','latex')
figure(1)
name1 = strcat('1d_dir_fd_csl_wn',wn,'_ppw',pts, ...
            '_pshift_',powershift,'_fshift_',factorshift,'.tex');
        
        if strcmp(pollution,'no')
            name1 = strcat('1d_dir_fd_csl_wn',wn,'_nopoll', ...
            '_pshift_',powershift,'_fshift_',factorshift,'.tex');
        end

matlab2tikz('filename',name1,'standalone',true); %save .tex file of tikz figure%

%End of eigenvalues CSL
close all;

%Eigenvalues of Def CSL
plot(real(eDCSL),imag(eDCSL),'.b','MarkerSize',16); 
axis equal
axis([-1 1.0 -1 1]);
xlabel('Re(z)','FontSize',14);
ylabel('Im(z)','FontSize',14);
set(gca,'Xtick',[-1 -0.5 0 0.5 1],'FontSize',14);
set(gca,'Ytick',[-1 -0.5 0 0.5 1],'FontSize',14);
hold on; 
hold on; 
t=linspace(0,2*pi,100);
plot(1/2+1/2*cos(t),1/2*sin(t),'k');

f = figure;

%setting file names for the eps images
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
name1 = strcat('1d_dir_fd_dcsl','_wn',wn,'_ppw',pts,...
            '_pshift_',powershift,'_fshift_',factorshift,'.tex');

        if strcmp(pollution,'no')
            name1 = strcat('1d_dir_fd_dcsl_wn',wn,'_nopoll', ...
            '_pshift_',powershift,'_fshift_',factorshift,'.tex');
        end
        
        
matlab2tikz('filename',name1,'standalone',true); %save .tex file of tikz figure%
end
