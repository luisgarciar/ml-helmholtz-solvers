function [err] = interfacefittedPoissonfemrate(pde,option,varargin)


%% Parameters
option = femoption(option); %#ok<*ASGLU>
maxIt = option.maxIt;
hmax = 0.01;

%% Initialize err
errL2 = zeros(maxIt,1);   errH1 = zeros(maxIt,1); 
erruIuh = zeros(maxIt,1); errMax = zeros(maxIt,1);
N = zeros(maxIt,1);

% generate a body-fitted mesh
[node,telem,selem,bdEdge] = interfacemesh(pde.phi,pde.box,hmax);


%% Finite Element Method        
for k = 1:maxIt

    [~,~,interfaceEdge] = findinterfaceedge(node,telem,pde.phi);
    % solve the equation
    [u,w,AE,AI,isExteriorTElem,isExteriorSEelm] = interfacefittedPoisson(node,telem,selem,pde,interfaceEdge,bdEdge, option);
    N(k) = size(node,1);
     
    NN = size(node,1);
    % compute error
    isInNode = false(NN,1);
    
    exteriorTElem = telem(isExteriorTElem,:);
    interiorTElem = telem(~isExteriorTElem,:);
    exteriorSElem = selem(isExteriorSEelm,:);
    interiorSElem = selem(~isExteriorSEelm,:);
    
    isInNode([interiorTElem(:);interiorSElem(:)]) = true;
    inNode = find(isInNode);
    u1 = u;
    u1(inNode) = u1(inNode) - w(inNode);
    
    e1 = getL2error(node,exteriorTElem,pde.exactuplus,u);
    e2 = getL2error(node,interiorTElem,pde.exactuminus,u1);
    e3 = getL2errorQ1(node,exteriorSElem,pde.exactuplus,u);
    e4 = getL2errorQ1(node,interiorSElem,pde.exactuminus,u1);
    errL2(k) = sqrt(e1^2+e2^2+e3^2+e4^2);
    
    e1 = getH1error(node,exteriorTElem,pde.Duplus, u,pde.d);
    e2 = getH1error(node,interiorTElem,pde.Duminus, u1,pde.d);
    e3 = getH1errorQ1(node,exteriorSElem,pde.Duplus, u,pde.d);
    e4 = getH1errorQ1(node,interiorSElem,pde.Duminus, u1,pde.d);
    errH1(k) = sqrt(e1^2+e2^2);
        
    uI = zeros(NN,1);
    uI1 = zeros(NN,1);
    idxE = [exteriorTElem(:);exteriorSElem(:)];
    idxI = [interiorTElem(:);interiorSElem(:)];
    uI(idxE) = pde.exactuplus(node(idxE,:)); % nodal interpolation
    uI1(idxI) = pde.exactuminus(node(idxI,:));
    erruIuh(k) = sqrt((u-uI)'*AE*(u-uI) + (u1-uI1)'*AI*(u1 - uI1));
    errMax(k) = max(max(abs(u(idxE)-uI(idxE)),max(abs(u1(idxI)-uI1(idxI)))));
    
    if k < maxIt
        [node,telem,selem,bdEdge] = interfacemesh(pde.phi,pde.box,hmax/2^k);
    end
    
end

%% Plot solution
if option.plotflag
     showsolution(node,exteriorTElem,u);
     hold on
    showsolution(node,exteriorSElem,u);
     showsolution(node,interiorTElem,u1);
     showsolution(node,interiorSElem,u1);
     hold off
end

%% Plot convergence rates
if option.rateflag
    figure;
    set(gcf,'Units','normal'); 
    set(gcf,'Position',[0.25,0.25,0.55,0.4]);
    subplot(1,2,1)
    showrate2(N(1:maxIt),errH1(1:maxIt),1,'-*','||Du-Du_h||',...
              N(1:maxIt),errL2(1:maxIt),1,'k-+','||u-u_h||');
    subplot(1,2,2)
    showrate2(N(1:maxIt),erruIuh(1:maxIt),1,'m-+','||Du_I-Du_h||',...
              N(1:maxIt),errMax(1:maxIt),1,'r-*','||u_I-u_h||_{\infty}');
end

%% Output
err = struct('N',N,'H1',errH1(1:maxIt),'L2',errL2(1:maxIt),...
             'uIuhH1',erruIuh(1:maxIt),'uIuhMax',errMax(1:maxIt));
