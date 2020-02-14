function [x,flag,relres,iter] = mlfgmres(b,x0,ml_mat,ml_prec,ml_prec_split,restrict,interp,maxiter)
%% MLFGMRES: Computes the solution of Ax=b by flexible GMRES via coarse
% grid iteration with multilevel deflation
%
% Input:
% ml_mat,ml_prec,ml_prec_split,restrict,interp
%%
numlevs = len(ml_mat);
n1 = size(ml_mat{1});
n2 = size(ml_prec{1});

assert(n1==n2,'inconsistent size of matrix and preconditioner');
n = n1;

%if numlevs==1 solve exactly and return
if numlevs == 1
    x = ml_mat{1}\b;
end

A = @(x) (ml_mat{1}*x);

npre  = 1; npos = 1; w = 0.5; smo = 'wjac';
x1    = zeros(length(A),1);
M     = @(x) feval(@Vcycle,ml_prec,ml_prec_split,restrict,interp,...
                    x1,x,npre,npos,w,smo,1);
maxit = maxiter(1);

%FGMRES initialization
r = b-A*x0; %initial residual
normr = norm(r);
normb = norm(b);

%Check tolerance of initial residual
if (normr < tol*normb)
    x      = x0;
    flag   = [];
    relres = [];
    iter   = [];
    return;
end

%Arnoldi basis V, matrix Z and upper Hessenberg H  s.t. AZ = VH
V = zeros(n,maxit+1);
Z = zeros(n,maxit+1);
H = zeros(maxit+1,maxit);

%right hand side of reduced system
g = norm(r)*eye(maxit+1,1);

%Givens rotation parameters
c = zeros(maxit+1,1);
s = zeros(maxit+1,1);

%initialization of Arnoldi basis
V(:,1) = r/norm(r);
j = 0;

while (j < maxit)
    j = j +1; 
    %Preconditioning step
    %the preconditioner is the multilevel deflation method
    %((Minv(I-A*P*(A_Hinv)*P^T) + P*inv(A_Hinv)*P^T))*V(:,j)
    z1 = restrict{1}*V(:,j);
    
    if numlevs == 2
        rH = ml_mat{2}\z1;
    else
        xinit = zeros(length(ml_mat{2}),1);
        %recursive call to multilevel method
        rH = mlfgmres(b,xinit,ml_mat{2:end},ml_prec{2:end},ml_prec_split,...
                      restrict{2:end},interp{2:end},maxiter{2:end});
    end
    
    eH = interp{1}*rH;
    z2 = V(:,j)-A*eH;
    z3 = feval(M,z2);
    Z(:,j) = z3 + eH;
    %end of preconditioning step
    
    
    w = A*Z(:,j);
    %Arnoldi Orthogonalization with Modified Gram-Schmidt
    for i=1:j
        H(i,j) = V(:,i)'*w;
        w = w - H(i,j)*V(:,i);
    end
    
    %Construction of new Arnoldi vector
    H(j+1,j) = norm(w);
    
    %Transform j-th column of Hessenberg matrix with previous Givens-rotations
    for i=1:j-1
        H(i:i+1,j) = [c(i),s(i); -s(i)',c(i)]*H(i:i+1,j);
    end
    
    beta = norm(H(j:j+1));   %Is this correct?
    if beta~=0
        V(:,j+1)=w/H(j+1,j); %Is this correct?
    %else %'happy breakdown %end
   
    %Compute new Givens rotation
        [s(j),c(j)]=givens(H(j,j),H(j+1,j));
    %Apply new Givens rotation to H
         H(j,j)   = cs(j)*H(j,j) + sn(j)*H(j+1,j);
         H(j+1,j) = 0.0;
    %Apply new Givens rotation to right hand side (Saad's book, eq 6.46)
         temp     =   c(j)*g(j);                   
         g(j+1)   =  -s(j)'*g(j);
         g(j)     =   temp;
    end
    
    error  = abs(g(j+1));
   
    if error <= tol 
        y = H(1:i,1:i) \ g(1:i); 
        x = x0 + Z(:,1:i)*y;
	    break 
    end
 
end
   

end


