function [g_n] = g_north(x,k,t)
%G_NORTH Function for testing the discretization
%of the Helmholtz equation with Sommerfeld 
%boundary conditions
%
% -div(grad u) - k^2u = 0 in Omega = (0,1)x(0,1)
%     du/dn - iku    = g on d(Omega). 
%
% g_n is g on the north boundary for the problem with
% exact solution
% u_ex(x,y) = cos(kxcost + kysint) + 1i*sin(kxcost + kysint)
%
%
%See Model problem 2 in 
%'Finite Element Analysis of Acoustic Scattering', Ihlenburg.
s   = sin(t); c = cos(t); 
xx  = k*(x*c+s);
g_n = k*(1-s)*(sin(xx)-1i*cos(xx));



