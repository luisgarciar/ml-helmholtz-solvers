%Test of 2D discretization of Poisson and Helmholtz problem
%We compare the discrete solution with a given exact solution
% of 
%
% -div(grad u) -k^2u = f in Omega = (0,1)x(0,1)
%                  u = 0 on boundary(Omega). 
%
% If u(x,y)=sin(m*pi*x)*sin(n*pi*y) then
% f =(-k^2+m^2*pi^2+pi^2*n^2)*sin(m*pi*x).*sin(n*pi*y)

clear all; clc; close all;
%k=0;  %Poisson problem
k = 10 %Helmholtz problem
m = 2; n=4; 
bc = 'dir';
u = @(x,y) sin(m*pi*x).*sin(n*pi*y);
f = @(x,y) (-k^2+m^2*pi^2+pi^2*n^2)*sin(m*pi*x).*sin(n*pi*y); %f=-u''-k^2u;

%% %% gridsize fixed  
% np = 25;
% npx = np; npy=np; hx = 1/(npx+1); hy = 1/(npy+1);
% xgrid = hx*(1:1:npx); ygrid = hy*(1:1:npy);
% np = npx*npy;
% 
% [X,Y]= meshgrid(xgrid,ygrid);
% u_ex = u(X,Y); u_ex = u_ex'; u_ex = reshape(u_ex,[np,1]); %exact solution
% b    = f(X,Y); b = b'; b = reshape(b,[np,1]);  %right hand side
% 
% A = helmholtz2(k,0,npx,npy,bc); b_d = A*u_ex; 
% u_d = A\b; 
% relerr1 = norm(u_d-u_ex,Inf)/norm(u_ex,Inf)
% 
% %Plotting the exact solution
% U_ex = reshape(u_ex,[npx,npy]);
% figure(1); surf(X,Y,U_ex); title('exact solution');
% 
% %Plotting the discrete rhs
% %b_d = reshape(b_d,[npx,npy]); 
% %figure(3); surf(X,Y,b_d); title('disc rhs')
% 
% %Plotting the exact rhs
% %b = reshape(b,[npx,npy]); 
% %figure(4); surf(X,Y,b); title('ex rhs')
% 
% 
% %Plotting the computed solution
% figure(2);
% U_d = reshape(u_d,[npx,npy]); 
% surf(X,Y,U_d); title('disc solution')
% 
% 
% %Plotting the error solution
% figure(3);
% Err = reshape(u_d-u_ex,[npx,npy]); 
% surf(X,Y,Err)


%% Refining gridsize, h->0

npt = 2.^(4:1:9);
relerr = zeros(length(npt),1);

for i=1:length(npt)
   nd = npt(i);
   npx = nd; npy=nd; hx = 1/(npx+1); hy = 1/(npy+1);
   xgrid = hx*(1:1:npx); ygrid = hy*(1:1:npy);
   np = npx*npy;
   [X,Y]= meshgrid(xgrid,ygrid);
   
   %exact solution and right hand side
   u_ex = u(X,Y); u_ex = u_ex'; u_ex = reshape(u_ex,[np,1]); %exact solution
   b    = f(X,Y); b = b'; b = reshape(b,[np,1]);  %right hand side
   
   %solving the discrete equation and computing error
   A = helmholtz2(k,0,npx,npy,bc); b_d = A*u_ex; 
   u_d = A\b; 
   relerr(i) = norm(u_d-u_ex,Inf)/norm(u_ex,Inf); 
end

loglog(npt,relerr)
title('relative error of solution vs number of points')
