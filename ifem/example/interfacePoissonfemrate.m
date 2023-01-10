function [err] = interfacePoissonfemrate(node,elem, pde,option,varargin)



%% Parameters
option = femoption(option); %#ok<*ASGLU>
maxIt = option.maxIt;

%% Initialize err
errL2 = zeros(maxIt,1);   errH1 = zeros(maxIt,1); 
erruIuh = zeros(maxIt,1); errMax = zeros(maxIt,1);
N = zeros(maxIt,1);

% generate a body-fitted mesh
[node,elem] = interfaceadaptivemesh(node,elem,pde.phi);

%% Finite Element Method        
for k = 1:maxIt

    [interiorElem,exteriorElem,interfaceEdge] = findinterfaceedge(node,elem,pde.phi);
    % solve the equation
    [u,w,AE,AI] = interfacePoisson(node,elem,pde,interfaceEdge,option);
    N(k) = size(node,1);
    
    NN = size(node,1);
    % compute error
    isInNode = false(NN,1);
    isInNode(interiorElem(:)) = true;
    inNode = find(isInNode);
    u1 = u;
    u1(inNode) = u1(inNode) - w(inNode);
    e1 = getL2error(node,exteriorElem,pde.exactuplus,u);
    e2 = getL2error(node,interiorElem,pde.exactuminus,u1);
    errL2(k) = sqrt(e1^2+e2^2);
    e1 = getH1error(node,exteriorElem,pde.Duplus, u);
    e2 = getH1error(node,interiorElem,pde.Duminus, u1);
    errH1(k) = sqrt(e1^2+e2^2);
        
    uI = zeros(NN,1);
    uI1 = zeros(NN,1);
    uI(exteriorElem(:)) = pde.exactuplus(node(exteriorElem(:),:)); % nodal interpolation
    uI1(interiorElem(:)) = pde.exactuminus(node(interiorElem(:),:));
    erruIuh(k) = sqrt((u-uI)'*AE*(u-uI) + (u1-uI1)'*AI*(u1 - uI1));
    errMax(k) = max(max(abs(u(exteriorElem(:))-uI(exteriorElem(:))),max(abs(u1(interiorElem(:))-uI1(interiorElem(:))))));
    
    if k < maxIt
        [node,elem] = interfaceUniformrefine(node,elem,interfaceEdge,pde.phi);
    end
    
end

%% Plot solution
if option.plotflag
    showsolution(node,exteriorElem,u);
    hold on
    showsolution(node,interiorElem,u1);
    hold off
end

%% Plot convergence rates
if option.rateflag
    figure;
    set(gcf,'Units','normal'); 
    set(gcf,'Position',[0.25,0.25,0.55,0.4]);
    subplot(1,2,1)
    showrate2(N(1:k),errH1(1:k),1,'-*','||Du-Du_h||',...
              N(1:k),errL2(1:k),1,'k-+','||u-u_h||');
    subplot(1,2,2)
    showrate2(N(1:k),erruIuh(1:k),1,'m-+','||Du_I-Du_h||',...
              N(1:k),errMax(1:k),1,'r-*','||u_I-u_h||_{\infty}');
end

%% Output
err = struct('N',N,'H1',errH1(1:k),'L2',errL2(1:k),...
             'uIuhH1',erruIuh(1:k),'uIuhMax',errMax(1:k));
