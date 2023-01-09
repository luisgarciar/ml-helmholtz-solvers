% Plots for Helmholtz book
clear all;
close all;

% Parameters
k   = 100;
ppw = 10;

b1  = 1; b2=0.5;
np = ceil(ppw*k/(2*pi))-1;
dim = 1;
bc = 'dir';
if (mod(np+1,2)==1) 
    np = np+1; 
end
npc = (np-1)/2;

%Eigenvalues of CSL and DCSL computed symbolically
[eigCSL,eigDCSL] = eigSL(k,ppw,b1,b2); 

%% % Eigenvalues of CSL,DCSL computed numerically
% A = helmholtz(k,np,bc);
% M = shift_laplace(k,b1,b2,np,bc);
% S = M\A;
% eigS = eig(full(S)); 
% 
% eigS = sort(eigS); %Eigenvalues of CSL
% 
% Y = fwrestriction(np,dim);
% Z = lininterpol(npc,dim);
% E = Z'*S*Z;   %coarse grid operator
% [L,U] = lu(E); I = eye(length(S));
% P = I-S*Z*(U\(L\Z'));
% PS = P*S;  %deflated operator;
% 
% eigPS = eig(PS); %Eigenvalues of DCSL computed numerically
% eigPS = sort(eigPS,'descend'); %Eigenvalues of DCSL computed numerically
% eigPS = eigPS(1:length(eigDCSL),1);

%% Plots to check things
%Checking that eigCSL and eigS coincide: They do coincide!
% figure(1); plot(real(eigCSL),imag(eigCSL),'+k');
% %hold on
% figure(2); plot(real(eigS),imag(eigS),'*r'); 

% %Checking that eigPA and eigDCSL coincide: They do coincide! (Test with small k)
% plot(real(eigPS),imag(eigPS),'*r'); 
% hold on
% plot(real(eigDCSL),imag(eigDCSL),'+k');
% axis equal

%% Plots for book
figure(1); plot(real(eigDCSL),imag(eigDCSL),'.b','MarkerSize',22);
axis square; %grid on

axis([0 1 -0.5 0.5])
xlabel('real(\lambda)','FontSize',14)
ylabel('imag(\lambda)','FontSize',14)

set(gca,'Xtick',[-1 -0.5 0 0.5  1],'FontSize',14);
set(gca,'Ytick',[-1 -0.5  0 0.5 1],'FontSize',14);

hold on;

% To construct the circle
angle = 0:0.005:2*pi;
r   = 0.5;    % Radius
xr  = 0.5; yr = 0;   % Center
xp=r*cos(angle);
yp=r*sin(angle);
plot(xr+xp,yr+yp,'-k')

figure(2); plot(real(eigCSL),imag(eigCSL),'.b','MarkerSize',22);
axis square; %grid on

axis([0 1 -0.5 0.5])

xlabel('real(\lambda)','FontSize',14)
ylabel('imag(\lambda)','FontSize',14)

set(gca,'Xtick',[-1 -0.5 0 0.5  1],'FontSize',14);
set(gca,'Ytick',[-1 -0.5  0 0.5 1],'FontSize',14);

hold on;

% To construct the circle
angle = 0:0.005:2*pi;
r   = 0.5;    % Radius
xr  = 0.5; yr = 0;   % Center
xp=r*cos(angle);
yp=r*sin(angle);
plot(xr+xp,yr+yp,'-k')

%% Old stuff
% axis square; grid on
% axis equal

% f = figure; plot(real(eDCSL),imag(eDCSL),'ok') 
% 
% % figure; plot(real(Lambda),imag(Lambda),'*k') 
% axis square; grid on
% axis equal
% axis([0 1.0 -0.5 0.5])
% xlabel('real(\lambda)','FontSize',14)
% ylabel('imag(\lambda)','FontSize',14)
% 
% set(gca,'Xtick',[0 0.2 0.4 0.6 0.8 1],'FontSize',14);
% set(gca,'Ytick',[-0.5 0 0.5],'FontSize',14);
% 
% hold on;
% 
% eigPA = eig(PA);
% eigA  = eig(A);
% 
% %tikzCSL(k,ppw,b1,b2)
%  

% f = figure; plot(real(eCSL),imag(eCSL),'ok')
% f = figure; plot(real(eDCSL),imag(eDCSL),'ok') 
% 
% % figure; plot(real(Lambda),imag(Lambda),'*k') 
% axis square; grid on
% axis equal
% axis([0 1.0 -0.5 0.5])
% xlabel('real(\lambda)','FontSize',14)
% ylabel('imag(\lambda)','FontSize',14)
% 
% set(gca,'Xtick',[0 0.2 0.4 0.6 0.8 1],'FontSize',14);
% set(gca,'Ytick',[-0.5 0 0.5],'FontSize',14);
% 
% hold on;
% 
% % To construct the circle
% angle = 0:0.01:2*pi;
% r   = 0.5;    % Radius
% xr  = 0.5; yr = 0;   % Center
% xp=r*cos(angle);
% yp=r*sin(angle);
% plot(xr+xp,yr+yp,'-k')
% 
% f.GraphicsSmoothing = 'off';
% ax = gca;           % get current axes
% ax.FontSmoothing = 'off'; 
% 
% 
% %  %setting file names for the eps images
% %  wn     = num2str(k);  pts = num2str(ppw);
% %  rshift = num2str(b1); ishift = num2str(b2);
% %  
% % %filename format: wavenumber_pointswavelength_realshift_imagshift.tex
% %  plot(real(eCSL),imag(eCSL),'b.');
% %  axis equal
% %  x = xlabel('$\mathrm{Re}(\lambda)$'); % x-axis label
% %  set(x,'Interpreter','latex')
% % 
% %  y=ylabel('$\mathrm{Im}(\lambda)$','interpreter','latex'); % x-axis label
% %  set(y,'Interpreter','latex')
% % 
% %  
% %  % figure(1)
% % % name1 = strcat('csl','_',wn,'_',pts,'_',rshift,'_',ishift,'.tex');
% % % matlab2tikz('filename',name1); %save .tex file of tikz figure%
% % % 
% % % plot(real(eDCSL),imag(eDCSL),'b.');
% % % axis equal
% % % xlabel('Re') % x-axis label
% % % ylabel('Im') % y-axis label
% % % set(gca, 'FontSize', 12)
% % %  
% % % name2 = strcat('dcsl','_',wn,'_',pts,'_',rshift,'_',ishift,'.tex');
% % % matlab2tikz(name2); %save .tex file of tikz figure%
% % 
% % %close(gcf)
