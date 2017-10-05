%Test of GMRES for the 2D Helmholtz problem with P1 finite elements
%
% We use the test problem
%
%  -div(grad u)-k^2 u   = f   in Omega= (0,1)x(0,1)
%  grad(u) dot n - i*ku = g on boundary(Omega)
%
% for the preconditioner we use the shifted Laplacian
%
%  -div(grad u)-(k^2 + i*eps) u   = f   in Omega= (0,1)x(0,1)
%  grad(u) dot n - i*ku = g on boundary(Omega)
%
%% Fixed wavenumber k and variable shift eps
kk   = [10 20 40];
itercsl = zeros(length(kk),1);
time_lu = zeros(length(kk),1);

for i=1:length(kk)
    k = kk(i);
    dim = 2;
    pollution = 'no';
    npf = ceil(k^(3/2));
    
    if (mod(npf+1,2)==0)  %set an odd number of interior points in 1D
        npf = npf+1;
    end
    npc = (npf-1)/2;
    
    poweps    = 1;
    factoreps = 1;
    bc = 'som';
    %Construct square mesh of meshsize h
    h = 1/npf;
    [node,elem] = squaremesh([0,1,0,1],h);
       
    %Find boundary nodes
    [bdNode,bdEdge,isBdNode] = findboundary(elem);
    
    %Sets Sommerfeld boundary conditions on all boundary edges
    bdFlag = setboundary(node,elem,'ABC');
    
    %The structures pde(helm,SL) contain data for the Helmholtz and 
    %shifted Laplace problems
    pdehelm = helmholtz2Dconstantwndata(k,0,1);
    pdeSL   = helmholtz2Dconstantwndata(k,factoreps,poweps);
    
    option.tol = 1e-12;
    [eqn1,~] = helmholtz2Dfem(node,elem,pdehelm,bdFlag,bdEdge);
    [eqn2,~] = helmholtz2Dfem(node,elem,pdeSL,bdFlag,bdEdge);
    
    %Helmholtz and shifted Laplace matrices
    A    = eqn1.A;
    Aeps = eqn2.A;
    
    %Parameters for GMRES
    restart   = [];
    tol       = 1e-8;
    maxit     = 200;
    
    %Setting up the preconditioner
    tic
    [L,U] = lu(Aeps);
    time_lu(i) = toc;
    AAepsinv  = @(x) A*(U\(L\x));  %Shifted Laplace preconditioned matrix
    b  = ones(length(A),1);
    x0 = zeros(size(b));
    
    %GMRES run
    [~, ~, ~, iter, ~] = gmres(AAepsinv, b, restart, tol, maxit);
    itercsl(i) = iter(2);
    
end