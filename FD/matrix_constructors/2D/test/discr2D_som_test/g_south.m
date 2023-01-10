function [g_s] = g_south(x,k,t)
%G_SOUTH Function for testing the discretization
% of the Helmholtz equation with Sommerfeld 
%boundary conditions
% -div(grad u) - k^2u = 0 in Omega = (0,1)x(0,1)
%     du/dn - iku    = g on d(Omega). 
%
% g_s is g on the south boundary for the problem with
% exact solution
% u_ex(x,y) = cos(kxcost + kysint) + 1i*sin(kxcost + kysint)
%
%
%See Model problem 2 in 
%'Finite Element Analysis of Acoustic Scattering', Ihlenburg.
%
s   = sin(t); c = cos(t);
g_s = k*(1+s)*(sin(k*x*c)-1i*cos(k*x*c));