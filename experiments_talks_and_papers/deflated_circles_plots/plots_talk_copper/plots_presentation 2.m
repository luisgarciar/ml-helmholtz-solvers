% Plots for presentation
k     = 10;       %wavenumber
pmin  = 7;       %min points per wavelenght
beta  = 0.5;      %complex shift
bc    = 'som' ; 
np    = ceil(k*pmin/(2*pi))-1; % heuristic for number of discretization points (TODO: find reference)
%kvar = @(x,y) k*ones(size(x));
%f    = @(x,y) 0*ones(size(x)); 

% %  Sommerfeld 1D Problem, spectrum k=10
%  [A,~,~] = helmholtz1d(f,k,np,bc,0);
%  [M,~,~] = shift_laplace1d(f,k,1,beta,np,bc,0);
%  
% S = M\A;
% eigvS = eig(full(S));
% eigvA = eig(full(A));
% 
% figure(2); plot(real(eigvA),imag(eigvA),'*b','MarkerSize',12);
% figure(1); plot(real(eigvS),imag(eigvS),'*b','MarkerSize',12);
% 
% hold on
% 
  angle = 0:0.01:2*pi;
r   = 0.5;    % Radius
%  xr  = 0.5; yr = 0;   % Center
%  xp=r*cos(angle);
%  yp=r*sin(angle);
%  plot(xr+xp,yr+yp,'-k')
 
% 
% %% Construction of Sommerfeld 2D matrices (See Ernst, Golub for notation)
 e = ones(np+2,1);
 h = 1/(np+1);
 
%% Helmholtz matrix
 T = spdiags([-e (4-k^2*h^2)*e -e], -1:1, np+2, np+2);
 T(1,1) = 3-k^2*h^2+1i*k*h; T(np+2,np+2) = 3-k^2*h^2+1i*k*h;
 a = -1+1i*k*h;
 
 J = spdiags([e 0*e e], -1:1, np+2, np+2);
 V = zeros(np+2,np+2); V(1,1)=1; V(np+2,np+2)=1;
 Id = speye(np+2);
 
  A = kron(Id,T)+kron(a*V-J,Id); %Helmholtz matrix
  A = A/h^2;
  eigvA = eig(full(A));
  
%Shifted Laplacian
 T      = spdiags([-e (4-k^2*(1-1i*beta)*h^2)*e -e], -1:1, np+2, np+2);
 T(1,1) = 3-k^2*(1-1i*beta)*h^2+1i*k*h; T(np+2,np+2) = 3-k^2*(1-1i*beta)*h^2+1i*k*h;
 a      = -1+1i*k*h; 
 J  = spdiags([e 0*e e], -1:1, np+2, np+2);
 V  = zeros(np+2,np+2); V(1,1)=1; V(np+2,np+2)=1;
 Id = speye(np+2);
 M  = kron(Id,T)+kron(a*V-J,Id); %Shifted Laplacian
 M = M/h^2;
 eigv = eig(full(M));

%figure(1); plot(real(eigv),imag(eigv),'*b','MarkerSize',12);%

S = M\A;
Sinv = A\M;
eigvS = eig(full(S));
figure(1); plot(real(eigvA),imag(eigvA),'.b','MarkerSize',20);
figure(1); plot(real(eigvS),imag(eigvS),'.b','MarkerSize',20);%
r = ceil(length(A)/2);
R = 10+5*rand(size(A));
[U N]=qr(R);
E = Z'*S*Z;  Id = speye(length(A));
P = Id-S*Z*(E\S');
PS = P*S;
eigvPS = eig(PS);

hold on

axis([-1 1 -1 1]);

set(gca,'Xtick',[-1 -0.5 0 0.5 1],'FontSize',14);
set(gca,'Ytick',[-1 -0.5 0 0.5 1],'FontSize',14);


angle = 0:0.01:2*pi;
r   = 0.5;    % Radius
xr  = 0.5; yr = 0;   % Center
xp  = r*cos(angle);
yp  = r*sin(angle);
plot(xr+xp,yr+yp,'-k');

xlabel('Re(\lambda)','FontSize',14)
ylabel('Im(\lambda)','FontSize',14)

matlab2tikz('filename','helmk30ppw12som.tex','standalone',true, 'extraaxisoptions',['xlabel style={font=\LARGE},'...
                       'ylabel style={font=\LARGE},','ticklabel style={font=\Huge}']);