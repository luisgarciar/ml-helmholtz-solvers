clear all;
close all;

k  = 1000;
n  = ceil(6*k/pi);
N  = 2*n+1;
h  = 1/(N+1);
b1 = 1;
b2 = 0.5;

%% Symbolic computation of the eigenvalues of Minv*P*A

cj  = cos((h*pi/2)*(1:n));
sj  = sin((h*pi/2)*(1:n));

g1 = (h*pi/2)*(1:n);
g2 = (h*pi/2)*(N+1-(1:n));

%Eigenvalues of Helmholtz operator
lj  = (4/h^2)*(sin(g1).^2)-k^2;
lNj = (4/h^2)*(sin(g2).^2)-k^2;

%Eigenvalues of shifted Laplacian
muj  = (4/h^2)*(sin(g1).^2) - k^2*(b1-1i*b2);
muNj = (4/h^2)*(sin(g2).^2) - k^2*(b1-1i*b2);

%Eigenvalues of preconditioned Shifted Laplacian
Slj  = lj./muj;
SlNj = lNj./muNj;

mj = lj.*(cj.^4)+lNj.*(sj.^4); 
aj = (cj.^4)./(sj.^4+cj.^4);

x       = (1./lNj).*aj + (1./lj).*(1-aj);
eigvPA  = (1./x);  % Eigenvalues of deflated Helmholtz
figure(1);
scatter(real(eigvPA),imag(eigvPA),'r+');

y  = (1./SlNj).*aj + (1./Slj).*(1-aj);
eigvTLMG = (1./y);
 
scatter(real(eigvTLMG),imag(eigvTLMG),'r+');
axis([0 2 -1 1 ])


%scatter(real(eigvPA),imag(eigvPA),'r+');

%scatter(real(mj),imag(mj),'r+');
%hold on



%% Checking the eigenvalues by constructing the matrices
% 
% 
%  A = helmholtz(k,N,'dir');
%  M = shift_laplace(k,b1,b2,N,'dir');
%  
%  P  = lininterpol(n);    %prolongation operator
%  R  = 0.5*P';            %restriction operator
%  E  = R*A*P;             %coarse grid matrix
%  D  = eye(N)-A*P*(E\R);  %deflation operator
%  Q  = P*(E\R);
%  D  = D + Q;
%  S  = (D*A);            %deflated matrix
%  %S  = M\(D*A);         %ADEF-1 preconditioned matrix
%  %D1 = D1+Q;
%  %S1 = D*(M\A);         %ADEF-1 preconditioned matrix
%  
% %eigvS = eig(full(S));
% %  eigvE = eig(full(E));
% 
% %figure(2);
% %scatter(real(eigvS),imag(eigvS),'b*');
% %  %scatter(real(eigvE),imag(eigvE),'b*');
% % 
% %  eigvA = eig(full(A));
% %  eigvsA=(4/h^2)*(sin(pi*h/2*(1:N)).^2) - k^2;
% %  eigvsA=eigvsA'
% %  
% %  scatter(real(eigvA),imag(eigvA),'r+');
% %  hold on
% %  scatter(real(eigvsA),imag(eigvsA),'b+');
% %  
% % %eigvA  = sort(eigvA);
% % %eigvsA = sort(eigvsA);
% %   
% %  %axis([0 2 -1 1 ])
% 
% %% We check the eigenvalues of the coarse grid system numerically 
% % and by constructing the matrices
% 
% %First we check the aliasing properties of the interpolation and fw 
% %restriction operators
% 
% %Fourier matrices of size N and n
% 
% %Fh= 


%lininterpol(5)
 
 
 
 
 
 
