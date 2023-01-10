function [p,u,eqn,info] = DarcyRT0(node,elem,pde,bdFlag,option)
%% DARCYRT0 Darcy equation: lowest order RT element.
%
%  [p,u] = DarcyRT0(node,elem,pde,bdFlag) produces an approximation of
%  the Poisson equation 
%
%       -div(d*grad(p))=f  in \Omega, with 
%       Dirichlet boundary condition p=g_D on \Gamma_D, 
%       Neumann boundary condition   d*grad(p)*n=g_N on \Gamma_N
%
%  The velocity u = d*grad(p) is approximated using the lowest order
%  Raviart-Thomas element and p by piecewise constant element.
%
%  [p,u] = PoissonRT0(node,elem,pde,bdFlag,option) specifies options
%   - option.solver
%     'direct': the built in direct solver \ (mldivide)
%     'tripremixPoisson':     multigrid-type solvers mg is used.
%     'uzawapcg': PCG for the Schur complement equation
%     'none': only assemble the matrix equation but not solve
%
%   The default setting is to use the direct solver for small size problems
%   and transforming based multigrid solvers for large size problems. 
%
%  It is almost identical to PoissonRT0 except (sigma, u) is changed to
%  (u,p) and the tensor is given by pde.K.
%
%  Example
%
%
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if ~exist('option','var'), option = []; end

time = cputime;  % record assembling time
%% Diffusion coefficient
if ~isfield(pde,'K'), pde.K = []; end
if isfield(pde,'K') && ~isempty(pde.K)
   if isnumeric(pde.K)
      K = pde.K;                   % K is an array
   else                            % K is a function
      center = (node(elem(:,1),:) + node(elem(:,2),:) + node(elem(:,3),:))/3;
      K = pde.K(center);  % take inverse sequencil.             
   end
else
    K = [];
end

%% Data structure
elemold = elem;
[elem,bdFlag] = sortelem(elem,bdFlag);  % ascend ordering
[elem2edge,edge] = dofedge(elem);
NT = size(elem,1); NE = size(edge,1);
[Dlambda,area,elemSign] = gradbasis(node,elem);

%% Assemble matrix 
Nu = NE; Np = NT; Ndof = Nu + Np;

% M. Mass matrix for RT0 element
M = getmassmatvec(elem2edge,area,Dlambda,'RT0',K);

% B. negative divergence operator
B = icdmat(double(elem2edge),elemSign*[1 -1 1]);

% C. zero matrix.
C = sparse(Np,Np);

A = [M B';B C];

%% Assemble right hand side.
fu = zeros(Np,1);
if ~isfield(pde,'f') || (isreal(pde.f) && (pde.f==0))
    pde.f = [];
end
if ~isfield(option,'fquadorder')
    option.fquadorder = 2;   % default order
end
if ~isempty(pde.f)
	[lambda,weight] = quadpts(option.fquadorder);
	nQuad = size(lambda,1);
    for k = 1:nQuad
		% quadrature points in the x-y coordinate
		pxy = lambda(k,1)*node(elem(:,1),:) ...
			+ lambda(k,2)*node(elem(:,2),:) ...
			+ lambda(k,3)*node(elem(:,3),:);
		fp = pde.f(pxy);
		fu = fu - fp*weight(k);
    end
    fu = fu.*area;
end
clear fp
F((Nu+1):Ndof,1) = fu;

%% Boundary Conditions
if ~exist('bdFlag','var'), bdFlag = []; end
[AD,F,bigu,freeDof,freeEdge,isPureNeumannBC] = getbdRT0(F);
eqn = struct('M',AD(1:NE,1:NE),'B',AD(NE+1:end,1:NE),'C',AD(NE+1:end,NE+1:end),...
             'f',F(1:NE),'g',F(NE+1:end),'freeDof',freeDof);

%% Record assembling time
assembleTime = cputime - time;
if ~isfield(option,'printlevel'), option.printlevel = 1; end
if option.printlevel >= 2
    fprintf('Time to assemble matrix equation %4.2g s\n',assembleTime);
end

%% Solve the linear system.
% Set up solver type
if isempty(option) || ~isfield(option,'solver')    % no option.solver
    if Ndof <= 2e3  % Direct solver for small size systems
        option.solver = 'direct';
    else         % MGCG  solver for large size systems
        option.solver = 'tripremixPoisson';
    end
% elseif strcmp(option.solver,'mg')
%     option.solver = 'tripremixPoisson';    
end
solver = option.solver;
% solve
switch lower(solver);
    case 'direct'
        t = cputime;
        bigu(freeDof) = AD(freeDof,freeDof)\F(freeDof);
        u = bigu(1:NE);
        p = bigu(NE+1:end); 
        info = struct('solverTime',cputime - t,'itStep',1,'error',0,'flag',0,'stopErr',0);
    case 'none'
        u = zeros(NE,1); p = zeros(NT,1); info =[];        
    case 'tripremixpoisson'
        [u,p,info] = tripremixPoisson(eqn.M,eqn.B,eqn.C,eqn.f,eqn.g,elemold);    
    case 'uzawapcg'
        [u,p,info] = uzawapcg(eqn.M,eqn.B,eqn.C,eqn.f,eqn.g,elemold);
    case 'mg'
        option.freeEdge = freeEdge;
        [u0,p,info] = mgDarcy(eqn.M,eqn.B,eqn.f,eqn.g,elemold,option);
        u = bigu(1:NE);
        u(freeEdge) = u0;
        Au = eqn.B*(eqn.B)';
        b = eqn.B*(eqn.f - eqn.M*u);
%         u = mg(Au,b,elemold);
        p = zeros(Np,1);
        if isPureNeumannBC
            freep = 1:Np-1;
        else
            freep = 1:Np;
        end
        p(freep) = Au(freep,freep)\b(freep);
end
if isPureNeumannBC % post process for u for pure Neumann boundary condition
    pbar = sum(p.*area)/sum(area);
    p = p - pbar;
end

%% Output information
info.assembleTime = assembleTime; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction getbdRT0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [AD,F,bigu,freeDof,freeEdge,isPureNeumannBC] = getbdRT0(F)
    %% GETBDRT0 Boundary conditions for Poisson equation: RT0 element.
    %
    %  Created by Ming Wang. Improved the check of edgeSign by Long Chen.

    bigu = zeros(Ndof,1);
    
    %% No boundary conditions
    if ~isfield(pde,'g_D'), pde.g_D = []; end
    if ~isfield(pde,'g_N'), pde.g_N = []; end

    %% Set up bdFlag
    if isempty(bdFlag) % no bdFlag information
       if ~isempty(pde.g_N) % case: Neumann
           bdFlag = setboundary(node,elem,'Neumann');
       elseif ~isempty(pde.g_D) % case: Dirichlet
           bdFlag = setboundary(node,elem,'Dirichlet');
       end
    end

    %% Find Dirichlet and Neumann dofs 
    if ~isempty(bdFlag)
        isDirichlet(elem2edge(bdFlag(:)==1)) = true;
        isNeumann(elem2edge(bdFlag(:)==2)) = true;
        % Direction of boundary edges may not be the outwards normal
        % direction of the domain. edgeSign is introduced to record this
        % inconsistency.
        edgeSign = ones(NE,1);
        idx = (bdFlag(:,1) ~= 0) & (elemSign == -1);% first edge is on boundary
        edgeSign(elem2edge(idx,1)) = -1;
        idx = (bdFlag(:,2) ~= 0) & (elemSign == 1); % second edge is on boundary
        edgeSign(elem2edge(idx,2)) = -1;
        idx = (bdFlag(:,3) ~= 0) & (elemSign == -1);% third edge is on boundary
        edgeSign(elem2edge(idx,3)) = -1;
    end
    Dirichlet = edge(isDirichlet,:);
    Neumann = edge(isNeumann,:); 
    isBdDof = false(Ndof,1); 
    isBdDof(isNeumann) = true;   % for mixed method, Neumann edges are fixed
    freeDof = find(~isBdDof);
    isFreeEdge = true(NE,1);
    isFreeEdge(isNeumann) = false;
    freeEdge = find(isFreeEdge);
    
    %% Dirichlet boundary condition (Neumann BC in mixed form)
    %   We need only modify the rhs on dof associated with Dirichlet
    %   boundary. Compute the int_e \Phi\cdot n g_D on the boundary using
    %   quadrature rules.
    if ~isempty(pde.g_D) && isnumeric(pde.g_D) && (pde.g_D==0)
        pde.g_D = [];
    end
    if ~isempty(pde.g_D) && any(isDirichlet) 
        if ~isfield(option,'gNquadorder')
            option.gNquadorder = 2;   % default order exact for linear gN
        end
        [lambda,weight] = quadpts1(option.gNquadorder);
        nQuad = size(lambda,1);
        for ip = 1:nQuad
        	pxy = lambda(ip,1)*node(Dirichlet(:,1),:)+...
                  lambda(ip,2)*node(Dirichlet(:,2),:);               
            F(isDirichlet) = F(isDirichlet) + weight(ip)*pde.g_D(pxy);
        end
        F(isDirichlet) = F(isDirichlet).*edgeSign(isDirichlet);
        % no edge length since the basis of sigma contains it.
    end

    %% Neumann boundary condition (Dirichlet BC in mixed form)
    if ~isempty(pde.g_N) && any(isNeumann)
        % modify the rhs to include Dirichlet boundary condition 
        mid = 1/2*(node(Neumann(:,1),:)+node(Neumann(:,2),:));
        ve = node(Neumann(:,1),:)-node(Neumann(:,2),:);
%         edgeLength = sqrt(sum(ve.^2,2)); 
        ne = [ve(:,2) -ve(:,1)]; % rotation of tangential vector
        if isnumeric(pde.g_N)
            evalg_N = pde.g_N;
        else
            evalg_N = pde.g_N(mid,ne);
        end
        bigu(isNeumann) = evalg_N.*edgeSign(isNeumann);
        F = F - A*bigu;
        F(isNeumann) = bigu(isNeumann);
    end
    
    %% Pure Neumann boundary condition
    isPureNeumannBC = false;
    if ~any(isDirichlet) && any(isNeumann)
        freeDof = freeDof(1:end-1);  % eliminate the kernel by enforcing u(NT) = 0;
        isBdDof(end) = true;
        isPureNeumannBC = true;
        F(NE+1:end) = F(NE+1:end) - mean(F(NE+1:end)); % normalize
%         F(end) = 0;
    end

    %% Modify the matrix
    %  Build Neumann boundary condition(Dirichlet BC in mixed form) into the
    %  matrix AD by enforcing  |AD(bdNode,bdNode)=I, 
    %  AD(bdNode,FreeNode)=0, AD(FreeNode,bdNode)=0|.
    if any(isBdDof)
       bdidx = zeros(Ndof,1); 
       bdidx(isBdDof) = 1;
       Tbd = spdiags(bdidx,0,Ndof,Ndof);
       T = spdiags(1-bdidx,0,Ndof,Ndof);
       AD = T*A*T + Tbd;
    else
       AD = A;
    end
    end
end