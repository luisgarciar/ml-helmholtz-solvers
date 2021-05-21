function [x,flag,resid,iter] = mlpfgmres(rhs,x0,ml_mat,ml_prec,ml_prec_split,restrict,interp,maxiter,tol)
%% MLPFGMRES: Computes the solution of Ax=b by flexible GMRES via coarse
% grid iteration with multilevel deflation
%
% Input:
% ml_mat,ml_prec,ml_prec_split,restrict,interp
%%
numlevs = length(ml_mat);
n1 = length(ml_mat{1});
n2 = length(ml_prec{1});

assert(n1==n2,'inconsistent size of matrix and preconditioner');
n = n1;
tol_exit=tol;
N = size(rhs,1);


%if numlevs==1 solve exactly and return
if numlevs == 1
    x = ml_mat{1}\rhs;
else
    
    A = @(x) (ml_mat{1}*x);
    
    %setup for shifted Laplace preconditioner
    x1    = zeros(length(rhs),1);
    npre  = 1; npos = 1; w = 2/3; smo = 'wjac';
    M     = @(x) feval(@Vcycle,ml_prec,ml_prec_split,restrict,interp,x1,x,npre,npos,w,smo,1);
    
    %FGMRES initialization
    maxit = maxiter(1);
    r = rhs-A(x0); %initial residual
    normr = norm(r);
    normb = norm(rhs);
    
    %Check tolerance of initial residual
    if (normr < tol*normb)
        x      = x0;
        flag   = [];
        resid  = [];
        iter   = [];
        return;
    end
    
    %Check if rhs == 0
    beta0 = normb;
    if (beta0==0)
        x = zeros(n,1);
        iter = 0;
        resid = 0;
        return;
    end
    

    %Initialize memory
    V = zeros(N, maxit+1);      %Krylov Subspace basis
    Z = zeros(N, maxit);        %Preconditioned basis
    H = zeros(maxit+1, maxit);  %Hessenberg matrix
    rhs = zeros(maxit+1,1);
    R = zeros(maxit,maxit);
    c = zeros(maxit,1);
    s = zeros(maxit,1);
    y = zeros(maxit,1);
    resid = zeros(maxit+1,1);
    
    %initialization of Krylov basis    
    if norm(x0)>0
        Ax0 = A(x0);
        r = rhs - Ax0;
    else
        r = rhs;
    end
    
    normr = norm(r);
    resid(1)=norm(r);
    V(:,1) = r/normr;
    b(1) = normr; %Right hand side for Hessenberg system
    i = 0;
    
    
    %Main iteration
    while (i < maxit)
        i = i +1;
        converge = 0;
       
        %Preconditioning step
        %the preconditioner is the multilevel deflation method
        %((Minv(I-A*P*(A_Hinv)*R) + P*inv(A_Hinv)*R))*V{j}
        z1 = restrict{1}*V(:,i);
        
        if numlevs == 2
            rH = ml_mat{2}\z1;
        else
            xinit = zeros(length(ml_mat{2}),1);
            %recursive call to multilevel method
            [rH,~,~,~] = mlpfgmres(z1,xinit,ml_mat(2:numlevs),ml_prec(2:numlevs),ml_prec_split(2:numlevs),restrict(2:numlevs-1),interp(2:numlevs-1),maxiter(2:end),tol);
        end
        
        eH = interp{1}*rH;
        z2 = v-A(eH);
        z3 = M(z2); % shifted Laplacian step
        Z(:,i) = z3 + eH;
        %end of preconditioning step
         
        %w = Az
        V(:,i+1) = A(Z(:,i));
        
        %Modified Gram-Schmidt
        for k = 1:i 
          H(k,i) = V(:,k)'*V(:,i+1);
          V(:,i+1) = V(:,i+1) - H(k,i)*V(:,k); 
        end % end for k
        
        % new orthonormal basis 
        H(i+1,i) = norm(V(:,i+1));
        V(:,i+1) = V(:,i+1)/H(i+1,i); % becareful small H(i+1,i)
        
         %--------------------------------------
        % Use Givens transformation to get upper triangular system R 
        %--------------------------------------
        R(1,i) = H(1,i);
        
        % apply the previous Givens transformations
        if (i~=1)
            for k = 2:i
                temp = c(k-1)*R(k-1,i) + s(k-1)*H(k,i);
                R(k,i) = -s(k-1)*R(k-1,i) + c(k-1)*H(k,i);
                R(k-1,i) = temp;           
            end          
        end 
        
        % new Givens transformation
        delta = sqrt(R(i,i)^2 + H(i+1,i)^2);
        c(i) = R(i,i)/delta;
        s(i) = H(i+1,i)/delta;
        
        R(i,i) = c(i)*R(i,i) + s(i)*H(i+1,i);
        
        % apply Givens transformation to Hessenberg right hand side b
        b(i+1) = -s(i)*b(i);
        b(i)   = c(i)*b(i);
        
        % count iterations
        iter  = i;
        
        % check convergence b(i+1) = || f-Au_k ||  
        resid(i) = abs(b(i+1));
        
        if ((resid(i)/resid(1)) < tol_exit)
            converge = 1;
            flag = 0;
            break;
        end
        
    end
    
    % solve the upper trangular matrix
    y(1:i) = R(1:i, 1:i)\rhs(1:i);
    
    % solution
    x = x0 + Z(:,1:i)*y(1:i);
   
end


end

