%Test of 2D discretization of Helmholtz problem with Sommerfeld
%boundary conditions
%
%
% We compare the discrete solution with a given exact solution
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
bc = 'som1';
t  = pi/2; c = cos(t); s = sin(t);
u  = @(x,y) cos(k*x*c + k*y*s) + 1i*sin(k*x*c + k*y*s);
npcc = 4;
ppw = 0.5;

[npf,numlev] = fd_npc_to_npf(npcc,k,ppw);


%% gridsize fixed
npint1D = npf;
npxint = npint1D; npyint = npint1D; hx = 1/(npxint+1); hy = 1/(npyint+1);
xgrid  = hx*(0:1:npxint+1); ygrid = hy*(0:1:npxint+1);
npgrid = (npxint+2)*(npyint+2);
[X,Y]  = meshgrid(xgrid,ygrid);
u_ex   = u(X,Y); u_ex = u_ex'; u_ex = reshape(u_ex,[npgrid,1]); %exact solution

%right hand side
g = zeros(npgrid,1);

%Fill the right hand side vector with boundary data
%Corner points
%South-West corner (0,0)
g(1,1) = g_south(0,k,t)/hy + g_west(0,k,t)/hx;
%South-East corner (1,0)
g(npxint+2,1) = g_south(1,k,t)/hy + g_east(0,k,t)/hx;
%North-West corner (0,1)
n      = (npxint+2)*(npyint+1)+1;
g(n,1) = g_north(0,k,t)/hy+g_west(1,k,t)/hx;
%North_East corner (1,1)
g(npgrid,1) = g_north(1,k,t)/hy+g_east(1,k,t)/hx;

%Noncorner boundary points

%West boundary: (0,j*hy)
for j = 1:npyint
    n = j*(npxint+2)+1; %changing from (0,j) to global index n
    g(n,1) = g_west(j*hy,k,t)/hx;
end

%South boundary: (0,j*hy)
for ii=1:npxint
    n=ii+1; %changing from (i,0) to global index n
    g(n,1) = g_south(ii*hx,k,t)/hy;
end

%East boundary: (1,j*hy)
for j=1:npyint
    n = (npxint+2)*(j+1); %changing from (nx,j) to global index
    g(n,1)=g_east(j*hy,k,t)/hx;
end

%North boundary: (i*hx,1)
for ii=1:npxint
    n = (npxint+2)*(npyint+1)+(ii+1);
    g(n,1) = g_north(ii*hx,k,t)/hy;
end
%Boundary data has been checked!

A = helmholtz2(k,0,npxint,npyint,bc);
%b_d = A*u_ex;
u_d = A\g;
relerr1 = norm(real(u_d)-real(u_ex),inf)/norm(real(u_ex),inf);

%Plotting the exact solution
U_ex = reshape(real(u_ex),[npxint+2,npyint+2]);
fig1 = figure(1); surf(X,Y,U_ex); title('exact solution');

%set(fig1, 'units', 'inches', 'position', [1 1 3 2])

% Plotting the computed solution
fig2 = figure(2);
U_d = reshape(real(u_d),[npxint+2,npyint+2]);
surf(X,Y,U_d); title('discrete solution')
%set(fig2, 'units', 'inches', 'position', [1 1 3 2])

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

length(relerr)


for ii=1:length(npt)
    nd = npt(ii);
    npxint = nd; npyint= nd; hx = 1/(npxint+1); hy = 1/(npyint+1);
    xgrid  = hx*(0:1:npxint+1); ygrid = hy*(0:1:npyint+1);
    [X,Y]  = meshgrid(xgrid,ygrid);
    npgrid = (npxint+2)*(npyint+2);
    
    %exact solution and right hand side
    u_ex = u(X,Y); u_ex = u_ex'; u_ex = reshape(u_ex,[npgrid,1]); %exact solution
    
    %Fill the right hand side vector with boundary data
    g = zeros(npgrid,1);
    %Corner points
    %South-West corner (0,0)
    g(1,1) = g_south(0,k,t)/hy + g_west(0,k,t)/hx;
    %South-East corner (1,0)
    g(npxint+2,1) = g_south(1,k,t)/hy + g_east(0,k,t)/hx;
    %North-West corner (0,1)
    n      = (npxint+2)*(npyint+1)+1;
    g(n,1) = g_north(0,k,t)/hy+g_west(1,k,t)/hx;
    %North_East corner (1,1)
    g(npgrid,1) = g_north(1,k,t)/hy+g_east(1,k,t)/hx;
    
    
    %Noncorner boundary points  
    %West boundary: (0,j*hy)
    for j = 1:npyint
        n = j*(npxint+2)+1; %changing from (0,j) to global index n
        g(n,1) = g_west(j*hy,k,t)/hx;
    end
    
    %South boundary: (0,j*hy)
    for i=1:npxint
        n=i+1; %changing from (i,0) to global index n
        g(n,1) = g_south(i*hx,k,t)/hy;
    end
    
    %East boundary: (1,j*hy)
    for j=1:npyint
        n = (npxint+2)*(j+1); %changing from (nx,j) to global index
        g(n,1)=g_east(j*hy,k,t)/hx;
    end
    
    %North boundary: (i*hx,1)
    for i=1:npxint
        n = (npxint+2)*(npyint+1)+(i+1);
        g(n,1) = g_north(i*hx,k,t)/hy;
    end
    
    %solving the discrete equation and computing error
    A = helmholtz2(k,0,npxint,npyint,bc); 
    u_d = A\g;
    relerr(ii) = norm(u_d-u_ex)/norm(u_ex);
    
end

loglog(npt,relerr)
title('relative error of solution vs number of points')
