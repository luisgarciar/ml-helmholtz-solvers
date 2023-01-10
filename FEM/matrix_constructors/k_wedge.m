function [k] = k_wedge(center,kref)
%K_WEDGE returns the space dependent wavenumber k(x,y)
%for the wedge problem in (0,1)x(0,1)
% Input: 
% (x,y): pair in (0,1)^2
%  kref: reference wavenumber (mid layer) 
%
% Output:
% k: space dependent wavenumber
%
% Author: Luis Garcia Ramos, 
%         Institut fur Mathematik, TU Berlin
%
%Version 0.1 - Nov 2015
%%%%%

x = center(1); y=center(2);

if(y>0 && y < 0.2*x+0.2)
    k=(4/3)*kref;
elseif(y >= 0.2*x+0.2 && y < -0.2*x+0.8)
    k = kref;
elseif(y >= -0.2*x+0.8 && y < 1)
    k = 2*kref;
else
    k = 0;
end
    
    
end

