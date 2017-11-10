function [eqn,info] = helmholtz2Dfem(node,elem,pde,bdFlag,bdEdge,option)
%% HELMHOLTZ2DFEM Helmholtz equation: P1 linear element.
%  Constructs a matrix for a 2D Helmholtz and Shifted Laplace model problem.
%  Constructs the matrix corresponding to the discretization of a 1D Helmholtz (or shifted Laplace)
%  problem with mixed or Sommerfeld boundary conditions
%
%   [eqn,info] = helmholtz2Dfem(node,elem,pde,bdFlag,option) produces the
%   matrix and right hand side corresponding  to the discretization of
%   the Helmholtz problem
%
%    -\div(grad u) - k^2u  = f in \Omega, with
%    -i*k*u + grad(u)*n    = g;  on \Gamma (Absorbing boundary condition)
%
%   or the shifted Laplace problem
%
%   -\div(grad u) - (k^2 + i*eps)u  = f in \Omega, with
%   -i*k*u + grad(u)*n = g  on \Gamma (ABC)
%
%   The mesh is given by node and elem and the info on boundary edges
%   is given by bdFlag. See meshdoc, bddoc for details.
%
%   The data is given by the
%   structure pde which contains function handles k, f, g
%   and the parameters poweps, factoreps for the shifted Laplacian
%   TO DO: Add parameters for the shifted Laplacian
%
%   The structure option contains the order of quadrature for the
%   integration of the right hand side function.
%
%   Input:
%
%   Example
%
%   Modified by Luis Garcia Ramos, based on the iFEM package
%   by Long Chen
%
%   TO DO: Fix problem with nonconstant wavenumbers in the
%   function getbd
%
%   Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if ~exist('option','var'), option = []; end
N = size(node,1); NT = size(elem,1);
dim = size(node,2);
Ndof = N;

if isempty(pde.poweps)
    factoreps=0;
else
    factoreps= pde.factoreps;
end

if isempty(pde.poweps)
    poweps=1;
else
    poweps= pde.poweps;
end

tic;  % record assembling time

%% Compute geometric quantities and gradient of local basis

%Dlambda(1:NT, 1:2, 1:3) is an array containing the gradients
%of local P1 basis functions
%For example, Dlambda(t,:,1) is the gradient of
% the basis function lambda corresponding to the 1st index (vertex) of
% the t-th triangle. The elemSign array taking values 1 or -1 records the
% sign of the signed area.

[Dphi,area] = gradbasis(node,elem);

%% Assemble stiffness matrix
% Delta is the negative Laplacian
% k2M is the squared wavenumber times the mass matrix
% M is the mass matrix
Delta   = sparse(Ndof,Ndof);
M       = sparse(Ndof,Ndof);
k2iepsM = sparse(Ndof,Ndof);

for i = 1:3
    for j = i:3
        % $A_{ij}|_{\tau} = \int_{\tau}\nabla \phi_i \cdot \nabla \phi_j dxdy$
        Aij = (Dphi(:,1,i).*Dphi(:,1,j) + Dphi(:,2,i).*Dphi(:,2,j)).*area;
        if (j==i)
            Delta = Delta + sparse(elem(:,i),elem(:,j),Aij,Ndof,Ndof);
        else
            Delta = Delta + sparse([elem(:,i);elem(:,j)],[elem(:,j);elem(:,i)],...
                [Aij; Aij],Ndof,Ndof);
        end
        %Assembling the mass matrix
        if isnumeric(pde.k)    %constant wave number
            k2     = (pde.k)^2;
            k2ieps = (pde.k)^2 + 1i*factoreps*(pde.k)^poweps;
        else %variable wave number k
            center = (node(elem(:,1),:) + node(elem(:,2),:) + node(elem(:,3),:))/3;
            k2     = (pde.k(center))^2;
            k2ieps = (pde.k(center))^2 + 1i*factoreps*pde.k(center)^poweps;
        end
        k2iepsMij = k2ieps*area*((i==j)+1)/12;
        Mij = area.*((i==j)+1)/12;
        if (j==i)
            k2iepsM = k2iepsM + sparse(elem(:,i),elem(:,j),k2iepsMij,Ndof,Ndof);
            M = M + + sparse(elem(:,i),elem(:,j),Mij,Ndof,Ndof);
        else
            k2iepsM = k2iepsM + sparse([elem(:,i);elem(:,j)],[elem(:,j);elem(:,i)],...
                [k2iepsMij; k2iepsMij],Ndof,Ndof);
            M =  M + sparse([elem(:,i);elem(:,j)],[elem(:,j);elem(:,i)],...
                [Mij; Mij],Ndof,Ndof);
        end
    end
end
clear Aij

A = Delta - k2iepsM;

%% Assemble the right hand side
b = zeros(Ndof,1);
if ~isfield(option,'fquadorder')
    option.fquadorder = 3; %default order of quadrature for integration
end
if ~isfield(pde,'f') || (isreal(pde.f) && (pde.f==0))
    pde.f = [];
end
if ~isempty(pde.f)
    [lambda,weight] = quadpts(option.fquadorder);
    phi = lambda;        % linear bases
    nQuad = size(lambda,1);
    bt = zeros(NT,3);
    for p = 1:nQuad
        % quadrature points in the x-y coordinate
        pxy = lambda(p,1)*node(elem(:,1),:) ...
            + lambda(p,2)*node(elem(:,2),:) ...
            + lambda(p,3)*node(elem(:,3),:);
        fp = pde.f(pxy); %evaluate f at the quadrature points
        for i = 1:3
            bt(:,i) = bt(:,i) + weight(p)*phi(p,i)*fp;
        end
    end
    bt = bt.*repmat(area,1,3);
    b = accumarray(elem(:),bt(:),[Ndof 1]);
end
clear pxy bt

%% Set up boundary conditions, see function getbd below
if ~exist('bdFlag','var'), bdFlag = []; end
[AD,b,freeNode] = getbd(b);

%% Record assembling time
assembleTime = toc;
if ~isfield(option,'printlevel'), option.printlevel = 1; end
if option.printlevel >= 2
    fprintf('Time to assemble matrix equation %4.2g s\n',assembleTime);
end

%% Output information
eqn = struct('A',AD,'b',b,'freeNode',freeNode,'Delta',Delta,'M',M);
info.assembleTime = assembleTime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions getbd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [AD,b,freeNode] = getbd(b)
        %% Set up of boundary conditions.
        
        % 1) Modify the right hand side b for the Sommerfeld boundary
        %   condition. The boundary integral of the function g is added to b.
        %
        % Reference: Long Chen. Finite Element Methods and its Programming.
        % Lecture Notes.
        
        %u = zeros(Ndof,1);
        %% Initial check
        if ~isfield(pde,'g'),   pde.g = 0; end
        if ~isfield(pde,'f'),   pde.f = 0; end
        
        %% Part 1: Modify the matrix for Sommerfeld BC
        % In the array bdFlag the Sommerfeld BC is marked as ABC:=3
        % This is a default option of iFEM, could be modified later.
        
        ABC = [];
        isABC = (bdFlag(:) == 3);
        if any(isABC)
            allEdge = [elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])];
            ABC = allEdge(isABC,:);
        end
        if ~isempty(ABC)
            ve = node(ABC(:,1),:) - node(ABC(:,2),:);
            edgeLength = sqrt(sum(ve.^2,2));
            %       mid = (node(ABC(:,1),:) + node(ABC(:,2),:))/2;
                        
            %Computation of boundary integral
            % int g phi_i phi_j ds
            ii = [ABC(:,1),ABC(:,1),ABC(:,2),ABC(:,2)];
            jj = [ABC(:,1),ABC(:,2),ABC(:,1),ABC(:,2)];
            temp = -sqrt(-1)*sqrt(k2)*edgeLength; 
            
            %computing the boundary integrals:
            %int_{edge} phi_i phi_j ds = 1/3*edgeLength if i==j,
            %int_{edge} phi_i phi_j ds = 1/6*edgeLength if i=/= j
            ss = [1/3*temp, 1/6*temp, 1/6*temp, 1/3*temp];
            A = A + sparse(ii,jj,ss,Ndof,Ndof); %?
             
            % Find Dirichlet boundary nodes: fixedNode
            fixedNode = []; freeNode = [];
            if ~isempty(bdFlag) % find boundary edges and boundary nodes
                [fixedNode,bdEdge,isBdNode] = findboundary(elem,bdFlag);
                freeNode = find(~isBdNode);
            end
            if isempty(bdFlag) && ~isempty(pde.g_D) && isempty(pde.g_N) && isempty(pde.g_R)
                % no bdFlag, only pde.g_D is given
                [fixedNode,bdEdge,isBdNode] = findboundary(elem);
                freeNode = find(~isBdNode);
            end
            
            if ~isempty(fixedNode)
                AD = A(freeNode,freeNode);
            else
                AD = A;
            end            
        end
        
        
        %% Part 2: Find boundary edges and modify the right hand side b
        
        %Old iFEM code: Neumann BCs (not considered here)
        % Find boundary edges: ABC
             if ~isempty(bdFlag)  % bdFlag specifies different bd conditions
                 ABC = bdEdge;  % bdEdge found in findboundary: line 275
             end
        %     if isempty(bdFlag) && ~isempty(pde.g_N)
        %         % no bdFlag, only pde.g_N or pde.g_R is given in the input
        %         [tempvar,Neumann] = findboundary(elem); %#ok<ASGLU>
        %     end
        
        % Absorbing boundary condition in the right hand side
        % In most cases, g= [] for Helmholtz equation
        if  isnumeric(pde.g) && all(pde.g == 0)
            pde.g = [];
        end
        
        if ~isempty(ABC) && ~isempty(pde.g)
            ve = node(ABC(:,1),:) - node(ABC(:,2),:);
            edgeLength = sqrt(sum(ve.^2,2));
            
            if ~isfield(option,'gquadorder')
                option.gquadorder = 2;   % default order exact for linear gN
            end
            [lambda_g,weight_g] = quadpts1(option.gquadorder);
            phi_g = lambda_g;                 % linear bases
            nQuadg = size(lambda_g,1);
            ge = zeros(size(ABC,1),2);
            for pp = 1:nQuadg
                % quadrature points in the x-y coordinate
                ppxy = lambda_g(pp,1)*node(ABC(:,1),:) ...
                    + lambda_g(pp,2)*node(ABC(:,2),:);
                gp = pde.g(ppxy);  %boundary function at the quadrature nodes
                for ig = 1:2
                    ge(:,ig) = ge(:,ig) + weight_g(pp)*phi_g(pp,ig)*gp;
                end
            end
            ge = ge.*repmat(edgeLength,1,2);%Jacobian (edgeLength)
            b = b + accumarray(ABC(:), ge(:),[Ndof,1]);
        end
        
        %     % The case with non-empty Neumann edges but g_N=0 or g_N=[] corresponds to
        %     % the zero flux boundary condition on Neumann edges and no modification of
        %     % A,u,b is needed.
        %
        %     % Dirichlet boundary condition
        %     if isnumeric(pde.g_D) && all(pde.g_D == 0)   % zero g_D
        %         pde.g_D = [];
        %     end
        %     if ~isempty(fixedNode) && ~isempty(pde.g_D)
        %         if isnumeric(pde.g_D)  % pde.g_D could be a numerical array
        %             u(fixedNode) = pde.g_D(fixedNode);
        %         else % pde.g_D is a function handle
        %             u(fixedNode) = pde.g_D(node(fixedNode,:));
        %         end
        %         b = b - A*u;
        % %         b(fixedNode) = u(fixedNode);
        %     end
        %     if ~isempty(fixedNode)
        %         b = b(freeNode);
        %     end
        %     % Difference with Poisson: no modification for pure Neumann case.
        %
        %     % The case with non-empty Dirichlet nodes but g_D=0 or g_D=[] corresponds
        %     % to the zero Dirichlet boundary condition and no modification of u,b is
        %     % needed.
    end % end of getbd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % end of Helmholtz