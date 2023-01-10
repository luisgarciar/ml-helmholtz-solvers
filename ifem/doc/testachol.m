% [node,elem] = squaremesh([0,1,0,1],1/16);
% load lakemesh
% load bunny
[node,elem] = cubemesh([0,1,0,1,0,1],0.125);
for k = 1:4
    [A,M] = assemblematrix(node,elem);
    A = A + 0.1*M;
%     A = A+200*M;
    b = ones(size(A,1),1);
    tol = 1e-6; maxit = 200;
    % tic; L1 = chol(A,'lower'); toc;
%     fprintf('\n Incomplete chol decomposition');
%     tic; 
%     L1 = ichol(A); 
%     R1 = L1';
%     toc;
%     fprintf('\n ichol as Preconditioner');
%     tic;
%     [x1,fl1,rr1,it1,rv1] = pcg(A,b,tol,maxit,@(x)icholpre(x,A,L1,R1));
%     toc;
%     fprintf('#dof: %8.0u,  iter: %2.0u\n',size(A,1), it1)
%     semilogy(0:it1,rv1./norm(b),'r.');
% 
    fprintf('\n Approximate chol decomposition');
    tic;
    [L2,p,Ac] = achol(A);
    toc;
    fprintf('\n Achol as Preconditioner');
    tic;
    B = tril(A);
    [x2,fl2,rr2,it2,rv2] = pcg(A,b,tol,maxit,@(r)acholpre(r,A,L2,L2',p,Ac,B,B',3));
    toc;
    fprintf('#dof: %8.0u,  iter: %2.0u\n',size(A,1), it2)
    
%     amg(A,b);
%     semilogy(0:it2,rv2./norm(b),'b.');
%     hold on

    [node,elem] = uniformrefine(node,elem);    
end