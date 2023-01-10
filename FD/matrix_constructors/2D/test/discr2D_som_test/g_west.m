function [g_w] = g_west(y,k,t)
%G_WEST Function for testing the discretization
%of the Helmholtz equation with Sommerfeld 
%boundary conditions
%
% -div(grad u) - k^2u = 0 in Omega = (0,1)x(0,1)
%     du/dn - iku     = g on d(Omega). 
%
% g_w is g on the west boundary for the problem with
% exact solution
% u_ex(x,y) = cos(kxcost + kysint) + 1i*sin(kxcost + kysint)
%
%See Model problem 2 in 
%'Finite Element Analysis of Acoustic Scattering', Ihlenburg.
%

s=sin(t); c=cos(t);
g_w = k*(c+1)*(sin(k*y*s)-1i*cos(k*y*s));

end

