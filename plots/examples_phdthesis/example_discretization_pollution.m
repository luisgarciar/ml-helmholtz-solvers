clc
clear global;
dim       = 2;
poweps    = 2;
factoreps = 0.5;

color1  = [0.5 0 0.5; 0 0 1; 0 0.5 0; 1 0 0];
kk = [10 20 40 60];


for s=1:length(kk)
    k    = kk(s);
    npc  = 3;
    par  = 0.9;
    npf  = ceil(k^(3/2));
    np   = npf-2;
    
    %% This section creates the data for the Helmholtz problem with exact solution
    %  -div(grad u)-k^2 u = 0 in Omega= (0,1)x(0,1)
    %  grad(u) dot n - i*ku = g in bd(Omega)
    %  with a plane wave as exact solution, i.e.,
    %  u = e^{i kk \cdot (x,y)}  where kk = k (cos(t), sin(t))
    %  and t is the angle (direction) of the wave
    %   Output:
    %   pde:    struct containing the following data:
    %   All function handles to be applied to input of size (N,2)
    %   'f':    function handle for right hand side (equal to 0)
    %   'exactu': function handle for exact solution
    %   'gradu': function handle for gradient of exact solution
    %   'k':  wavenumber
    %   'g': function handle for boundary data
    
    t = pi/2; % direction of propagation
    pdehelm = helmholtz2Dplanewavedata(k,t);
    u_ex    = pdehelm.exactu;
    Du = pdehelm.gradu;
    f = pdehelm.f;
    g = pdehelm.g;
    
    % Data for shifted Laplace problem (no need for right hand side)
    pdeSL      = helmholtz2Dconstantwndata(k,factoreps,poweps);
    option.tol = 1e-8;
    
    %% Constructing coarse grid
    [npf,numlev] = fem_npc_to_npf(npc,k,par);  %number of points in finest grid (1D)
    h = 1/(npc+1);
    [node,elem] = squaremesh([0 1 0 1],h);
    
    %L2error
    error_L2    = zeros(numlev+1,1);
    L2norm_uex  = zeros(numlev+1,1);
    relerror_L2 = zeros(numlev+1,1);
    
    %H1 error
    error_H1    = zeros(numlev+1,1);
    H1norm_uex  = zeros(numlev+1,1);
    relerror_H1 = zeros(numlev+1,1);
    
    
    %Refining the grid
    for j=1:numlev+1
        
        [node,elem] = uniformrefine(node,elem);
        [bdNode,bdEdge,isBdNode] = findboundary(elem);
        bdFlag = setboundary(node,elem,'ABC');
        
        % Assembly of the matrix and rhs vector
        % stored in struct eqn1
        [eqn1,~] = helmholtz2Dfem(node,elem,pdehelm,bdFlag,bdEdge);
        
        %Setup of multilevel Krylov solver and shifted Laplacian
        [mg_matHelm,~,~,~]= mg_setupfem_2D(npc,j+1,pdehelm);
        [mg_matCSL,mg_splitCSL,restrCSL,interpCSL]= mg_setupfem_2D(npc,j+1,pdeSL);
        
        %Helmholtz Matrix and shifted Laplace matrix
        A      = mg_matHelm{1};
        Aeps   = mg_matCSL{1};
        
        %Matrix and Right hand side
        b  = eqn1.b;
        
        if (length(b)<=100000)
            tic
            u = A\b;
            time_lu = toc;
        else
            
            %% Solving the Galerkin problem with the shifted Laplacian & multilevel FGMRES
            npre = 1; npos = 1; w  = 0.4; smo = 'wjac';
            u0      = zeros(length(A),1);
            %            Aepsinv = @(x) Vcycle(mg_matCSL,mg_splitCSL,restrCSL,interpCSL,u0,x,npre,npos,w,smo,1);
            %            AP  = @(x) A*Aepsinv(x);
            tol = 1e-8;
            restart = [];
            maxiter = ones(length(mg_matHelm),1);
            maxiter(1:5)=[200,8,4,2,2]';
            
            %Solving the problem with mlfgmres
            tic
            [u,~,~,iter3] = mlfgmres(b,u0,mg_matHelm,mg_matCSL,mg_splitCSL,restrCSL,interpCSL,maxiter,tol);
            time_mlk = toc;
            
        end
        
        %plotting the real part of the solution
        % figure(1)
        % showsolution(node,elem,real(full(u)));
        % reu_ex = @(x) real(u_ex(x));
        z = zeros(size(u));
        
        %L2error
        error_L2(j)    = getL2error(node,elem,u_ex,u);
        L2norm_uex(j)  = getL2error(node,elem,u_ex,z);
        relerror_L2(j) = error_L2(j)/L2norm_uex(j);
        
        reDu = @(x) real(Du(x));
        %H1 error
        error_H1(j)    = getH1error(node,elem,Du,u);
        H1norm_uex(j)  = getH1error(node,elem,Du,z);
        relerror_H1(j) = error_H1(j)/H1norm_uex(j);
        
    end
    
    hh  = (1/sqrt(2)*2.^(-2:-1:(-2-(numlev))))';
    label = strcat('$k = ', num2str(kk(s)),'$');
    plt = loglog(hh,relerror_L2,'Color',color1(s,:),'LineWidth',2,'DisplayName',label,'Marker','+'...
    ,'MarkerSize',5);
    hold on
   
    
end

L = legend('$k$ = 10','$k$ = 20','$k$ = 40','$k$ = 60');
set(L,'Interpreter','latex','FontSize',17,'Location','Northwest');
set(gca,'FontSize',14)



%Hkerror
% error_Hk   = sqrt(error_H1.^2+k^2.*error_L2^2);
% Hknorm_uex = sqrt(H1norm_uex^2+k^2*error_L2^2);
% relerrorHk = error_Hk/Hknorm_uex
%error = norm(u_iter3-u_gal,Inf)
%time_lu
%time_mlk
% Computing the relative error of the Galerkin problem w.r.t.
% the exact solution
% We need the
%plotting the real part of the solution
%figure(1)
%showsolution(node,elem,real(full(u)));




