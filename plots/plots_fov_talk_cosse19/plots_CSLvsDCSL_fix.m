clear all;
close all;

ppw = 30;  %points per wavelength
bc = 'som';

%Parameters for the shifted Laplacian
%We use a shifted Laplacian of the form
% M =  A-i*eps*I, where A is the Helmholtz matrix
% and eps = factoreps*k^poweps
%The standard choice is factoreps = 0.5, poweps=2


%Choose which plots are generated
plot_csl    = 'yes';
plot_defcsl = 'yes';

% Parameters
%kk = [40 80 100 150];
k = 200;
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
M = mass(npf,bc);
X = sqrtm(full(M));

Ahat = Aeps\A;
Ahats = X\(Ahat*X);
[fvCSL, eCSL] = fv(full(Ahats),1,32);

plot(real(fvCSL),imag(fvCSL),'k','MarkerSize',10,'LineWidth',2);
hold on
plot(real(eCSL),imag(eCSL),'r+','MarkerSize',5,'LineWidth',2);


axis equal
axis([0  1.4 -0.7 0.7]);
xlabel('Re(z)','FontSize',14);
ylabel('Im(z)','FontSize',14);
set(gca,'Xtick',[-1 -0.5 0 0.5 1],'FontSize',14);
set(gca,'Ytick',[-1 -0.5 0 0.5 1],'FontSize',14);
hold on;

%Tikz Axis formatting
x = xlabel('$\mathrm{Re}(z)$');
% x-axis label
set(x,'Interpreter','latex')
y=ylabel('$\mathrm{Im}(z)$','interpreter','latex'); % x-axis label
set(y,'Interpreter','latex')
figure(1)


%%
R = fwrestrictionfem(npf,dim,bc);
Z = R';     %Prolongation. Deflation subspace: columns of Z
dim_def = size(Z,2);
%Field of values of restricted matrix
%To compute the field of values of the restricted matrix
%(DEFCSL) we need a basis Y for the orthogonal complement of
%the columns of Z, i.e., the nullspace of R.
N = length(A);
[eCSL,eDCSL] = eigSL(k,N,eps);
plot(real(eDCSL),imag(eDCSL),'+b','MarkerSize',15);
axis([0  1.4 -0.7 0.7]);
xlabel('Re(z)','FontSize',14);
ylabel('Im(z)','FontSize',14);
set(gca,'Xtick',[-1 -0.5 0 0.5 1],'FontSize',14);
set(gca,'Ytick',[-1 -0.5 0 0.5 1],'FontSize',14);
axis equal

