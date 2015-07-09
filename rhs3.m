function f = rhs3(x,y)
% RHS3 Computes a test function for 
%      the Helmholtz equation
% INPUT: 
% [x,y]: meshgrid on (0,1)^2
% OUTPUT:
% f: eigenfunction evaluated on meshgrid x,y 
%%
k  = 0; 
f  = (-4-k^2)*ones(size(x));
end
