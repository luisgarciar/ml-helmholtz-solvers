function [x] = smoother(U,L,D,P,b,x0,w,numit,smo)
%% SMOOTHER Basic relaxation methods 
%  Applies numit iterations of a basic relaxation method
%  (Gauss-Seidel, w-Jacobi, Red-Black Gauss Seidel)
%  to solve Ax=b with initial guess x_0
%
%   Use: [x] = smoother(U,L,D,P,b,x0,w,numit,smo)
%
%   Input:
%
%   U, L, D: (strictly) upper,lower and diagonal parts of A
%         P: Permutation matrix for Red-Black Gauss-Seidel (not working!!)
%     x0, b: Initial guess, right hand side    
%     numit: number of iterations 
%         w: relaxation parameter for Jacobi iteration
%       smo: type of smoother 
%            'gs': Gauss Seidel
%            'wjac': w-Jacobi
%            'rbgs': red-black Gauss Seidel (for 2D problems - currently not working!!)
%
%  Author: Luis Garcia Ramos, 
%          Institut fur Mathematik, TU Berlin
%          Version 0.1, Jun 2016
%
% TO DO: Fix Red-Black Gauss-Seidel
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
        for i=1:numit            
            x = (D+L)\(-U*x+b);
        end   
        
    case  'wjac'
        for i=1:numit
            x = w*(D\(-L*x-U*x+b))+(1-w)*x;
            %x = w*(D\((-L-U)*x)+D\b)+(1-w)*x;
            
            %Dinv=(1./D);
            %x = x + w*D\(b-(D*x);
        end
        
    case  'rbgs'
        x = P*x0;  Pb=P*b;  %Permuting data
        PLP = P*L*P'; PUP=P*U*P'; PDP=P*D*P';
        M = PDP+PLP;  N = -PUP;
        
        for i=1:numit
            x = M\(N*x+Pb);
        end
        x = P'*x;
end

end