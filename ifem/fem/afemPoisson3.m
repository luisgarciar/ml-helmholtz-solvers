function [err,time,solver,eqn,node,elem] = afemPoisson3(node,elem,pde,bdFlag,option,HB,varargin)

%% Check input arguments
if nargin >=1 && ischar(node)
    option.elemType = node;
    clear node
end
if ~exist('node','var') || ~exist('elem','var')
    % default mesh: Lshape
    [node,elem] = cubemesh([-1,1,-1,1,-1,1],1);
    [node,elem] = delmesh(node,elem,'x>0 & y<0');
end
if ~exist('option','var'), option = []; end
if ~exist('pde','var')
    pde = Lshapedata3;                          % default data
end
if ~exist('bdFlag','var')
    bdFlag = setboundary3(node,elem,'Dirichlet');
end
if ~exist('HB','var')
    N0 = size(node,1);
    HB(1:N0,1:3) = repmat((1:N0)',1,3); 
    HB(1:N0,4) = 0;
end

%% Parameters
option = afemoption(option,3);
maxIt = option.maxIt;
refType = option.refType;
elemType = option.elemType;
theta = option.theta;

%% Initialize err
errL2 = zeros(maxIt,1);   errH1 = zeros(maxIt,1); erreta = zeros(maxIt,1);
erruIuh = zeros(maxIt,1); errMax = zeros(maxIt,1);
errTime = zeros(maxIt,1); solverTime = zeros(maxIt,1); 
assembleTime = zeros(maxIt,1); meshTime = zeros(maxIt,1); 
itStep = zeros(maxIt,1);  stopErr = zeros(maxIt,1); flag = zeros(maxIt,1);
N = zeros(maxIt,1);

%% Generate an initial mesh 
for k = 1:option.L0
    if strcmp(refType,'red')
        [node,elem,bdFlag,HB] = uniformrefine3(node,elem,bdFlag,HB);
    elseif strcmp(refType,'bisect')
        [node,elem,bdFlag,HB] = uniformbisect3(node,elem,bdFlag,HB);
    end
end

%%  Adaptive Finite Element Method
% *SOLVE* -> *ESTIMATE* -> *MARK* -> *REFINE*
for k = 1:maxIt
    % Step 1: SOLVE
    switch elemType
        case 'P1'     % piecewise linear function P1 element
            [u,Du,eqn,info] = Poisson3(node,elem,pde,bdFlag,option,HB);
        case 'CR'     % piecewise linear function CR element
            [u,Du,eqn,info] = Poisson3CR(node,elem,pde,bdFlag,option,HB);
        case 'P2'     % piecewise quadratic function
            [u,Du,eqn,info] = Poisson3P2(node,elem,pde,bdFlag,option,HB);
        case 'WG'     % weak Galerkin element
            [u,Du,eqn,info] = Poisson3WG(node,elem,pde,bdFlag,option,HB);            
    end
    % compute error
    tic;
    if isfield(pde,'Du') 
        if ~isempty(Du)
            errH1(k) = getH1error3(node,elem,pde.Du,Du);
        else
            errH1(k) = getH1error3(node,elem,pde.Du,u);            
        end
    end
    if isfield(pde,'exactu')
        errL2(k) = getL2error3(node,elem,pde.exactu,u);        
        % interpolation
        switch elemType
            case 'P1'
                uI = Lagrangeinterpolate(pde.exactu,node,elem);
            case 'CR'
                uI = Lagrangeinterpolate(pde.exactu,node,elem,'CR',eqn.face);
            case 'P2'
                uI = Lagrangeinterpolate(pde.exactu,node,elem,'P2',eqn.edge);
            case 'WG'
                uI = Lagrangeinterpolate(pde.exactu,node,elem,'WG',eqn.face);
        end
        erruIuh(k) = sqrt((u-uI)'*eqn.A*(u-uI));
        errMax(k) = max(abs(u-uI));
    end
    errTime(k) = toc;
    % record time
    solverTime(k) = info.solverTime;
    assembleTime(k) = info.assembleTime;
    if option.printlevel>1
        fprintf('Time to compute the error %4.2g s \n H1 err %4.2g    L2err %4.2g \n',...
            errTime(k),errH1(k), errL2(k));    
    end
    % record solver information
    itStep(k) = info.itStep;
    stopErr(k) = info.stopErr;
    flag(k) = info.flag;
    % plot 
    N(k) = size(node,1);
    if option.plotflag && N(k) < 1e4 % show mesh and solution for small size
       figure(1);  showresult3(node,elem,u,[],[74 20]);   
    end
    % Step 2: ESTIMATE
    switch option.estType
        case 'recovery' % recovery type
            eta = estimaterecovery3(node,elem,u);         
        case 'residual' % residual type
            switch elemType
                case 'P1'
                    eta = estimateresidual3(node,elem,u,pde,bdFlag);
                case 'WG'
                    eta = estimateresidual3WG(node,elem,Du,pde);                    
            end
    end
    erreta(k) = sqrt(sum(eta.^2));    
    % Step 3: MARK
    switch option.markType
        case 'L2'
            markedElem = mark(elem,eta,theta);
        case 'MAX'
            markedElem = mark(elem,eta,theta,'MAX');            
    end
    % Step 4: REFINE
    if N(k) > option.maxN
        break;
    end
    [node,elem,bdFlag,HB] = bisect3(node,elem,markedElem,bdFlag,HB);
end

%% Plot convergence rates
if option.rateflag
    figure;
    set(gcf,'Units','normal'); 
    set(gcf,'Position',[0.25,0.25,0.55,0.4]);
    subplot(1,2,1)
    showrate2(N(1:k),errH1(1:k),10,'-*','||Du-Du_h||',...
              N(1:k),errL2(1:k),10,'k-+','||u-u_h||');
    subplot(1,2,2)
    showrate2(N(1:k),erruIuh(1:k),10,'m-+','||Du_I-Du_h||',...
              N(1:k),errMax(1:k),10,'r-*','||u_I-u_h||_{\infty}');
end

%% Output
err = struct('N',N(1:k),'H1',errH1(1:k),'L2',errL2(1:k),...
             'uIuhH1',erruIuh(1:k),'uIuhMax',errMax(1:k),'eta',erreta(1:k));
time = struct('N',N(1:k),'err',errTime(1:k),'solver',solverTime(1:k), ...
              'assemble',assembleTime(1:k),'mesh',meshTime(1:k));
solver = struct('N',N(1:k),'itStep',itStep(1:k),'time',solverTime(1:k),...
                'stopErr',stopErr(1:k),'flag',flag(1:k));
            
%% Display error and CPU time
display('Table: Error')
colname = {'#Dof','||u-u_h||','||Du-Du_h||','||DuI-Du_h||','||uI-u_h||_{max}','eta'};
disptable(colname,err.N,[],err.L2,'%0.5e',err.H1,'%0.5e',err.uIuhH1,'%0.5e',...
                  err.uIuhMax,'%0.5e',err.eta,'%0.5e');
     
display('Table: CPU time')
colname = {'#Dof','Assemble','Solve','Error','Mesh'};
disptable(colname,time.N,[],time.assemble,'%0.2e',time.solver,'%0.2e',...
                  time.err,'%0.2e',time.mesh,'%0.2e');     