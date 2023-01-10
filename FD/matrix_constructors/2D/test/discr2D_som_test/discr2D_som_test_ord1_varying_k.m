%Test of 2D discretization of Helmholtz problem with Sommerfeld
%boundary conditions
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
bc = 'som1';
t  = pi/2; c = cos(t); s = sin(t);
npcc = 4;
ppw = 0.5;

%% Varying wavenumber, k->infty

wn = [10 20 40 60 80];
relerr = zeros(length(wn),1);

for ii=1:length(wn)
    k = wn(ii);
    t  = pi/2; c = cos(t); s = sin(t);
    u  = @(x,y) cos(k*x*c + k*y*s) + 1i*sin(k*x*c + k*y*s);
    [nd,~] = fd_npc_to_npf(npcc,k,ppw);
    
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
    u_h = A\g;
    relerr(ii) = norm(real(u_h)-real(u_ex),inf)/norm(real(u_ex),inf);
    
end

semilogy(wn,relerr)
title('relative error of solution vs wavenumber')
