clear all;
close all;

ppw = 20;  %points per wavelength
bc = 'som';

%Parameters for the shifted Laplacian
%We use a shifted Laplacian of the form
% M =  A-i*eps*I, where A is the Helmholtz matrix
% and eps = factoreps*k^poweps
%The standard choice is factoreps = 0.5, poweps=2

%Choose which plots are generated
plot_csl    = 'yes';

% Parameters
%kk = [40 80 100 150];
k = 40;
poweps     =  2;
factoreps  =  5;

close all
eps = factoreps*(k^poweps);

%if pollution = 'no', np = ceil(k^(3/2))
pollution = 'yes';
%otherwise the number of grid points is chosen
%with a fixed number of points per wavelength
npf = ceil(ppw*k/(2*pi))-1;
if strcmp(pollution,'no')
    npf = ceil(k^(3/2));
end

if (mod(npf,2)==0)
    npf = npf+1;
end
npc = (npf-1)/2;

dim = 1;

A = helmholtzfem(k,npf,0,bc); %Helmholtz matrix
Aeps = helmholtzfem(k,npf,eps,bc); %Shifted Laplace matrix

Ahat = Aeps\A;
eCSL = eig(full(Ahat));

n = 100;
t = 2*pi*linspace(0,1,n);
x = 0.5 + 0.5*cos(t);
y = 0.5*sin(t);

plot(x,y,'k','LineWidth',1);
hold on
plot(real(eCSL),imag(eCSL),'b.','MarkerSize',18,'LineWidth',2);


axis equal
axis([-0.2  1.2 -0.7 0.7]);
xlabel('Re(z)','FontSize',20);
ylabel('Im(z)','FontSize',20);
set(gca,'Xtick',[-1 -0.5 0 0.5 1],'FontSize',20);
set(gca,'Ytick',[-1 -0.5 0 0.5 1],'FontSize',20);
hold on;

%Tikz Axis formatting
x = xlabel('$\mathrm{Re}(z)$');
% x-axis label
set(x,'Interpreter','latex')
y=ylabel('$\mathrm{Im}(z)$','interpreter','latex'); % x-axis label
set(y,'Interpreter','latex')
figure(1)