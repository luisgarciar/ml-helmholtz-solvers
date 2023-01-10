function [g_e] = g_east(y,k,t)
%G_EAST Function for testing the discretization
% of the Helmholtz equation with Sommerfeld 
%boundary conditions
%
% -div(grad u) - k^2u = 0 in Omega = (0,1)x(0,1)
%     du/dn - iku    = g on d(Omega). 
%
% g_e is g on the east boundary for the problem with
% exact solution (plane wave in direction t)
% u_ex(x,y) = cos(kxcost + kysint) + 1i*sin(kxcost + kysint)%
%
% See Model problem 2 in 
%'Finite Element Analysis of Acoustic Scattering', Ihlenburg.
%
s  = sin(t); c = cos(t);
yy   = k*(y*s+c);
g_e  = k*(1-c)*(sin(yy)-1i*cos(yy));
