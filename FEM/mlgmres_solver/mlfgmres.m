function [x,flag,resvec,iter] = mlfgmres(b,x0,ml_mat,ml_prec,ml_prec_split,restrict,interp,maxiter,tol)
%% MLFGMRES: Computes the solution of Ax=b by flexible GMRES via coarse
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

%if numlevs==1 solve exactly and return
if numlevs == 1
    x = ml_mat{1}\b;
else
    
    A = @(x) (ml_mat{1}*x);
    
    %setup for shifted Laplace preconditioner
    x1    = zeros(length(b),1);
    npre  = 1; npos = 1; w = 2/3; smo = 'wjac';
    M     = @(x) feval(@Vcycle,ml_prec,ml_prec_split,restrict,interp,x1,x,npre,npos,w,smo,1);
    maxit = maxiter(1);
    
    %FGMRES initialization
    r = b-A(x0); %initial residual
    normr = norm(r);
    normb = norm(b);
    
    %Check tolerance of initial residual
    if (normr < tol*normb)
        x      = x0;
        flag   = [];
        resvec = [];
        iter   = [];
        return;
    end
    
    beta0 = norm(b);
    if (beta0==0)
        x = zeros(n,1);
        iter = 0;
        resvec = 0;
        return;
    end
    
    %resvec = zeros(maxit+1,1);
    % Householder data
    W = zeros(maxit+1,1);
    R = zeros(maxit, maxit);
    J = zeros(2,maxit);
    
    
    %Krylov Subspace basis
    V = cell(maxit, 1);
    % Preconditioned subspace
    Z = cell(maxit, 1);
    
    %initialization of Krylov basis
    
    if norm(x0)>0
        Ax0 = A(x0);
        r = b - Ax0;
    else
        r = b;
    end
    u = r;
    normr = norm(r);
    
    % Prepare the first Householder vector
    beta = scalarsign(r(1))*normr;
    u(1) = u(1)+beta;
    u = u/norm(u);
    V{1} = u;
    
    j = 0;
    resvec(1)=norm(r);
    %First Hessenberg entry: Apply Householder projection to r
    W(1) = -beta;
    
    while (j < maxit)
        j = j +1;
        
        %Form P1*P2*P3...Pj*ej.
        % v = Pj*ej = ej - 2*u*u'*ej
        v = -2*u(j)'*u;
        v(j) = v(j) + 1;
        % v = P1*P2*...Pj-1*(Pj*ej)
        for i = (j-1):-1:1
            v = v - 2*V{i}*(V{i}'*v);
        end
        % Explicitly normalize v to reduce the effects of round-off.
        v = v/norm(v);
        
        %Preconditioning step
        %the preconditioner is the multilevel deflation method
        %((Minv(I-A*P*(A_Hinv)*R) + P*inv(A_Hinv)*R))*V{j}
        z1 = restrict{1}*v;
        
        if numlevs == 2
            rH = ml_mat{2}\z1;
        else
            xinit = zeros(length(ml_mat{2}),1);
            %recursive call to multilevel method
            [rH,~,~,~] = mlfgmres(z1,xinit,ml_mat(2:numlevs),ml_prec(2:numlevs),ml_prec_split(2:numlevs),restrict(2:numlevs-1),interp(2:numlevs-1),maxiter(2:end),tol);
        end
        
        eH = interp{1}*rH;
        z2 = v-A(eH);
        z3 = M(z2);
        Z{j} = z3 + eH;
        %end of preconditioning step
        w = A(Z{j});
        
        % Orthogonalize the Krylov vector
        %  Form Pj*Pj-1*...P1*Av.
        for i = 1:j
            w = w - 2*V{i}*(V{i}'*w);
        end
        
        % Update the rotators: Determine Pj+1.
        if (j ~= length(w))
            %  Construct u for Householder reflector Pj+1.
            u = [zeros(j,1); w(j+1:end)];
            alpha = norm(u);
            if (alpha ~= 0)
                alpha  = scalarsign(w(j+1))*alpha;
                u(j+1) = u(j+1)+alpha;
                u = u/norm(u);
                V{j+1} = u;
                
                %  Apply Pj+1 to v.
                %  v = v - 2*u*(u'*v);
                w(j+2:end) = 0;
                w(j+1) = -alpha;
            end
        end
        
        %  Apply Given's rotations to the newly formed v.
        for colJ = 1:j-1
            tmpv = w(colJ);
            w(colJ)   = conj(J(1,colJ))*w(colJ) + conj(J(2,colJ))*w(colJ+1);
            w(colJ+1) = -J(2,colJ)*tmpv + J(1,colJ)*w(colJ+1);
        end
        
        %  Compute Givens rotation Jm.
        if (j~=length(w))
            rho = norm(w(j:j+1));
            J(:,j) = w(j:j+1)./rho;
            W(j+1) = -J(2,j).*W(j);
            W(j) = conj(J(1,j)).*W(j);
            w(j) = rho;
            w(j+1) = 0;
        end
        
        
        R(:,j) = w(1:maxit);
        
        %Residual
        resid = abs(W(j+1))/beta0;
        resvec=[resvec', abs(W(j+1))]';
        
        if (resid<tol_exit)
            iter = j;
            flag = 0;
            break;
        end
    end
    
    % Correction
    y = R(1:j,1:j) \ W(1:j);
    dx = zeros(n,1);
    for i=j:-1:1
        dx = dx + y(i)*Z{i};
    end
    
    x = x0+dx;
    flag = 0;
    iter = j;
end


end


function sgn = scalarsign(d)
sgn = sign(d);
if (sgn == 0)
    sgn = 1;
end
end


