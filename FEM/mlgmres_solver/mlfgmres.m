function [x,r,Res,lmax,Z0, E0,count] = gmresdef(A,b,x0, tol,omg,grid,maxgrid, Prec0, lmax, maxiter,Z0, E0,count)
%Function to calculate the solution of Ax=b by GMRES-iteration via coarse
%grid iteration with deflation 

%actualization of the grid level, matrix size and input of iterations on
%the coarse level

[n,n]= size(A); % only square matrices allowed in the method -> change?
grid = grid+1;
iter = 0;

% Initialize some matrices and vectors
V  = sparse(zeros(n,maxiter(grid)+1));          % Arnoldi vectors
Z  = sparse(zeros(n,maxiter(grid)));
h  = sparse(zeros(maxiter(grid)+1,1));  	  % upper Hessenberg st A*V = V*H ...
QT = sparse(zeros(maxiter(grid)+1,maxiter(grid)+1));  % orthogonal factor st QT*H = R
R  = sparse(zeros(maxiter(grid),maxiter(grid)));	  % upper triangular factor st H = Q*R
f  = sparse(zeros(maxiter(grid),1));		  % y = R\f => x = x0 + V*y

r  = b-A*x0;

h(1)    = norm(r);
V(:,1)  = r/h(1);
QT(1,1) = 1;
phibar  = h(1);
normp0  = phibar;
error   = abs(phibar)/phibar;

%fprintf('%5.0f%18.8e%18.8e\n',iter,abs(phibar),error);
%size of Prec has to correspond to number of deflation vectors!

Prec=speye(floor(n/2));
% calculate the matrices only once, the first time they are needed and
% store them in cell arrays

if lmax(grid+1)==-1 
  % subspace deflation vector construction, depending only on matrix size
  %  Z0{grid}=subdef(n);

    % coarse grid matrix and next level preconditioner
    %%E0{grid}= Z0{grid}'*A*inv(Prec0)*Z0{grid}; 
    
   %  Prec=diag(diag(E0{grid}));
    
    % estimation for largest eigenvalue for later shift

   %lmax(grid+1)= omg*lmaxest(E0{grid}*inv(Prec)); 
   lmax(grid+1)= omg;

end

[sizeE0,sizeE0]=size(E0{grid});
% actual method

while iter < maxiter(grid) && error > tol
%
  iter = iter + 1;
%
  %Z(:,iter) = M\V(:,iter);
  z6 = A*(Prec0\V(:,iter));
  z5 = Z0{grid}'*(z6-lmax(grid)*V(:,iter));
  if grid < maxgrid
    [z4,~,~,lmax,Z0, E0,count] = gmresdef(E0{grid}, z5, zeros(sizeE0,1),tol, omg, grid, maxgrid, Prec, lmax, maxiter,Z0, E0,count);
  else
    z4 = E0{grid}\z5;
    count=count+1;
    %z4=gmresdirect(E0, zeros(sizeE0,1), z5, eye(sizeE0), [], 50, tol );
  end
  
  z3 = Z0{grid}*z4;
  Z(:,iter) = V(:,iter) - z3;
  z2 = A*(Prec0\z3);
  u = z6 - z2;
%
% Arnoldi orthogonalization
%  
  for k = 1:iter
    h(k) = V(:,k)'*u;
    u    = u - h(k)*V(:,k);
  end
%
% Construct the new Arnoldi vector
%
  h(iter+1)      = norm(u);
  V(:,iter+1)    = u/h(iter+1);
  R(1:iter,iter) = QT(1:iter,1:iter)*h(1:iter);
  rt             = R(iter,iter);
%
% Givens rotation
% 
  if h(iter+1) == 0
    c = 1.0;
    s = 0.0;
  elseif abs(h(iter+1)) > abs(rt)
    temp = rt/h(iter+1);
    s    = 1.0/sqrt(1.0 + abs(temp)^2);
    c    = -temp*s;
  else
    temp = h(iter+1)/rt;
    c    = 1.0/sqrt(1.0 + abs(temp)^2);
    s    = -temp*c;
  end
  
  R(iter,iter) = c'*rt - s'*h(iter+1);
  
  q = QT(iter,1:iter);
  QT(iter  ,1:iter) = c'*q;
  QT(iter+1,1:iter) = s*q;
  QT(iter  ,iter+1) = -s';
  QT(iter+1,iter+1) = c;
  f(iter)           = c'*phibar;
  phibar            = full(s*phibar);
  res(iter)         = abs(phibar);
%
  error = abs(phibar)/normp0;
  
  if grid==1 
       fprintf('%5.0f%18.8e%18.8e\n',iter,abs(phibar),error)
      % issparse(error)
 end
   
 % if (error <= tol)
 %   y = R(1:iter,1:iter)\f(1:iter);
 %   x = x0 + Prec0\(Z(:,1:iter)*y);
 % end
%    resnorm(iter,1) = norm(res_act);
%    resprec(iter,1) = norm(inv(M)*(b-A*x));
end
y = R(1:iter,1:iter)\f(1:iter);
x = x0 + Prec0\(Z(:,1:iter)*y);

if grid==1
    r=b-A*x;
end
%
%
Res=res/normp0;
    


end
