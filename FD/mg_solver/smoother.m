function [x0] = smoother(U,L,D,P,b,x0,w,numit,smo)
%% SMOOTHER Basic relaxation methods 
%  Applies numit iterations of a basic relaxation method
%  (Gauss-Seidel, w-Jacobi, Red-Black Gauss Seidel)
%  to solve Ax=b with initial guess x_0
%
%   Use: [x] = smoother_som(U,L,D,P,b,x0,w,numit,smo)
%
%   Input:
%   U, L, D: (strictly) upper,lower and diagonal parts of A=U+D+L
%         P: Permutation matrix for Red-Black Gauss Seidel
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

switch smo
    case 'gs'
        %M = D+L; N = -U;
        for i=1:numit            
            %x0 = M\(N*x0+b);
            x0 = (D+L)\(-U*x0+b);
        end   
        
    case  'wjac'
        for i=1:numit
            N = -(L+U); 
            x0 = w*(D\(N*x0+b))+(1-w)*x0;
            %x0 = w*(D\(-L*x0-U*x0+b))+(1-w)*x0;          
        end
        
    case  'rbgs'
        Mrb = P*(D+L)*P'; Nrb=-P*U*P';       
        x0 = P*x0; brb=P*b;        
        for i=1:numit            
            x0 = Mrb\(Nrb*x0+brb);
        end   
        x0 = P'*x0;
end

end