function [x] = smoother_som(U,L,D,P,b,x0,w,numit,smo)
%% SMOOTHER_SOM Basic relaxation methods 
%  Applies numit iterations of a basic relaxation method
%  (Gauss-Seidel, w-Jacobi, Red-Black Gauss Seidel)
%  to solve Ax=b with initial guess x_0
%
%   Use: [x] = smoother_som(U,L,D,P,b,x0,w,numit,smo)
%
%   Input:
%   U, L, D: (strictly) upper,lower and diagonal parts of A=U+D+L
%     x0, b: Initial guess, right hand side    
%     numit: number of iterations 
%         w: relaxation parameter for Jacobi iteration
%       smo: type of smoother 
%            'gs': Gauss-Seidel
%            'wjac': w-Jacobi
%            'rbgs': red-black Gauss-Seidel (for 2D problems - currently not working)
%
%  Author: Luis Garcia Ramos, 
%          Institut fur Mathematik, TU Berlin
%          Version 1.0, Sep 2016
%
%Includes Red-Black Gauss-Seidel
%

%%

[n,m]=size(U);
assert(n==m, 'Incorrect matrix size');
[n,m]=size(D);
assert(n==m, 'Incorrect matrix size');
[n,m]=size(L);
assert(n==m, 'Incorrect matrix size');

x = x0;
switch smo
    case 'gs'
        M = D+L; N = -U;
        for i=1:numit            
            x = M\(N*x+b);
        end   
        
    case  'wjac'
        for i=1:numit
            N = -(L+U); 
            x = w*(D\(N*x)+D\b)+(1-w)*x;
            %Dinv=(1./D);
            %x = x + w*D\(b-A*x);
        end
        
    case  'rbgs'
        Mrb = P*(D+L)*P'; Nrb=-P*U*P';       
        x = P*x; brb=P*b;        
        for i=1:numit            
            x = Mrb\(Nrb*x+brb);
        end   
        x = P'*x;

end

end