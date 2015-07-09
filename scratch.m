
clc
clear all

np = 10;
h  = 1/(np+1);
k  = 10;
dim = 1;

d       = -2*ones(np,1); 
l       =  ones(np,1)*(1); 
Lapl1d  =  spdiags([l d l],[-1 0 1],np,np); %Discrete Laplace operator

Ad_1    = -(1/h^2)*Lapl1d + k^2*speye(np);   %1D Discrete Helmholtz Operator


% if dim==1
%  x = h:h:1-h;
%  b = exp(x).*(x.*x+3*x+1);
%  u = x.*(1-x).*exp(x);
% end

x = h:h:1-h;
b = 2.*ones(size(x))';
u = (x.*(1-x))';

sol = Ad_1\b;
norm(sol-u);

error = zeros(7,1);
cond  = zeros(7,1);

np = ceil(10.^linspace(1,4,20));


for i=1:length(n)
    np = np(i);
    h  = 1/(np+1)
    
    % Construction of the matrix
    d       = -2*ones(np,1); 
    l       =  ones(np,1)*(1); 
    Lapl1d  =  spdiags([l d l],[-1 0 1],np,np); %Discrete Laplace operator
    Ad_1    = -(1/h^2)*Lapl1d- k^2*speye(np);   %1D Discrete Helmholtz Operator

    %Construction of the right hand side and the exact solution   
    x = h:h:1-h;
    b = 2.*ones(size(x))';
    u = (x.*(1-x))';

    sol = Ad_1\b;
    error(n,1) = norm(sol-u);
    cond(n,1)  = condest(Ad_1);
    
 
end


plot(log(error))

% 
% 
% %Dirichlet 2D matrix (only interior points)
% nv     = np^2;
% Lapl2d = kron(Lapl1d, speye(np)) + kron(speye(np),Lapl1d);
% Ad_2   =  -(1/h^2)*(Lapl2d- k^2*h^2*speye(nv));      %2D Discrete Helmholtz Operator
% 
% 


