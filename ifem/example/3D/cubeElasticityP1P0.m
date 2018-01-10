close all; 
% clear all;

%% Parameters 
maxIt = 3; 
N = zeros(maxIt,1);
errL2 = zeros(maxIt,1); 
% errH1 = zeros(maxIt,1);  
% erruIuh = zeros(maxIt,1);

para.mu = 0.35;
para.lambda = 0.2;
% pde = elasticitydata3(para);
pde = elasticity3datapoly;

%% Mesh
cube = [0 1 0 1 0 1];
[node,elem] = cubemesh(cube,1/2);
bdFlag = setboundary3(node,elem,'Dirichlet');

%% FEM
for k = 1:maxIt
%     bdFlag = setboundary3(node,elem,'Dirichlet','abs(z)>eps','Neumann','abs(z)<eps');
%     [uh,Du,eqn,info] = Poisson3Q1(node,elem,pde,bdFlag,option);
    [sigma,u,eqn,info] = elasticity3mfemP1P0(node,elem,pde,bdFlag);
    N(k) = size(elem,1);
    center = (node(elem(:,1),:) + node(elem(:,2),:) ...
            + node(elem(:,3),:) + node(elem(:,4),:))/4;
    uI = pde.exactu(center);
    e = uI(:) - u;    
    errL2(k) = max(abs(e));
%     errL2(k) = sqrt(e'*e/N(k));
%     errL2(k) = getL2error3(node,elem,pde.exactu,u);        
    if k < maxIt
        [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
    end
end

%% Plot convergence rates
figure;
showrate(N,errL2,2,'k-+','||u-u_h||');      
colname = {'#Dof','||u-u_h||'};
disptable(colname,N,[],errL2,'%0.5e');