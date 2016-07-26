%test of 1D discretization
k = 0; np = 30; h = 1/(np+1); a=2;   bc='dir'
u = @(x)  sin(a*pi*x);
f = @(x) ((a*pi)^2-k^2)*sin(a*pi*x); %f=-u''-k^2u;

npp = logspace(1,5,5); 
err =zeros(length(npp),1);
relerr = err;

for i=1:length(npp)
    np   = npp(i); h = 1/(np+1);
    grid = h*(1:1:np)'; b = f(grid);
    u_ex = u(grid);
    [A]  = helmholtz(k,np,bc);
    u_d  = A\b; 
    err(i) = norm(u_d-u_ex,Inf);
    relerr(i) = norm(u_d-u_ex,Inf)/norm(u_ex,Inf);
    plot(grid,u_d,'r');
    hold on
    plot(grid,u_ex,'b');
    pause(1)
    close all
end