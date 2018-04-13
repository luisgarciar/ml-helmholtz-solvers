%Test of 2D discretization of Poisson and Helmholtz problem
%We compare the discrete solution with a given exact solution
% of 
%
% -div(grad u) - k^2u = 0 in Omega = (0,1)x(0,1)
%  du/dn - iku        = g on d(Omega). 
%
% For t fixed, we consider the function (exact solution)
% u_ex(x,y)=cos(kxcost + kysint) + 1i*sin(kxcost + kysint)
%
clear all; clc; close all;
k  = 10; %Helmholtz problem
bc = 'som';
t  = pi/2; c = cos(t); s=sin(t);
u  = @(x,y)cos(k*x*c + k*y*s) + 1i*sin(k*x*c + k*y*s);


%% gridsize fixed  
 npint1D = 1000; 
 npxint = npint1D; npyint=npint1D; hx = 1/(npxint+1); hy = 1/(npyint+1);
 xgrid  = hx*(0:1:npxint+1); ygrid = hy*(0:1:npxint+1);
 npgrid = (npxint+2)*(npyint+2);
 [X,Y]  = meshgrid(xgrid,ygrid);
 u_ex   = u(X,Y); u_ex = u_ex'; u_ex = reshape(u_ex,[npgrid,1]); %exact solution

 %right hand side
 g = zeros(npgrid,1);
 
 %Fill the right hand side vector with boundary data
 %Corner points 
 %South-West corner (0,0)
 g(1,1) = 2*g_south(0,k,t)/hy+2*g_west(0,k,t)/hx;
 %South-East corner (1,0)
 g(npxint+2,1) = 2*g_south(1,k,t)/hy+2*g_east(0,k,t)/hx; 
 %North-West corner (0,1)
 n = (npxint+2)*(npyint+1)+1;
 g(n,1) = 2*g_north(0,k,t)/hy+2*g_west(1,k,t)/hx;
 %North_East corner (1,1)
 g(npgrid,1) = 2*g_north(1,k,t)/hy+2*g_east(1,k,t)/hx;
 
%Noncorner boundary points
%West boundary: (0,j*hy)
for j = 1:npyint
    n = j*(npxint+2)+1; %changing from (0,j) to global index n
    g(n,1) = 2*g_west(j*hy,k,t)/hx;          
end
        
%South boundary: (0,j*hy)
for i=1:npxint
    n=i+1; %changing from (i,0) to global index n
    g(n,1) = 2*g_south(i*hx,k,t)/hy; 
end

%East boundary: (1,j*hy)
for j=1:npyint
    n = (npxint+2)*(j+1); %changing from (nx,j) to global index
    g(n,1)=2*g_east(j*hy,k,t)/hx;    
end

%North boundary: (i*hx,1)
for i=1:npxint
    n = (npxint+2)*(npyint+1)+(i+1);
    g(n,1)=2*g_north(i*hx,k,t)/hy;
end
%Boundary data has been checked! 

A = helmholtz2(k,0,npxint,npyint,bc);
%b_d = A*u_ex; 
u_d = A\g; 
relerr1 = norm(real(u_d)-real(u_ex),inf)/norm(real(u_ex),inf);

%Plotting the exact solution
U_ex = reshape(real(u_ex),[npxint+2,npyint+2]);
figure(1); surf(X,Y,U_ex); title('exact solution');

% Plotting the computed solution
figure(2);
U_d = reshape(real(u_d),[npxint+2,npyint+2]); 
surf(X,Y,U_d); title('discrete solution')

% 
% %Plotting the discrete rhs
% %b_d = reshape(b_d,[npx,npy]); 
% %figure(3); surf(X,Y,b_d); title('disc rhs')
% 
% %Plotting the exact rhs
% %b = reshape(b,[npx,npy]); 
% %figure(4); surf(X,Y,b); title('ex rhs')
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
