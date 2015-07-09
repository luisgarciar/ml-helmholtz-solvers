function exsol = solrhs3(x,y)
% exsol Computes the exact solution for rhs3
% INPUT: 
% [x,y]: meshgrid on (0,1)^2
% OUTPUT:
% ex: exact solution for test problem
%%
exsol = x.*y.*(x-1).*(y-1); 
end