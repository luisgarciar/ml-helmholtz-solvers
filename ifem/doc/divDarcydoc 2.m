%% u = DIVDARCY
%
% Find a u such at Bu = g
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.


%% Set up 
% [node,elem] = squaremesh([0,1,0,1],0.125); 
[node,elem] = circlemesh(0,0,1,0.25);
[elem2edge,edge] = dofedge(elem);
showmesh(node,elem);

pde = sincosNeumanndata;
bdFlag = setboundary(node,elem,'Neumann');
option.solver = 'none';
[~,~,eqn] = PoissonRT0(node,elem,pde,bdFlag,option);
freeEdge = eqn.freeEdge;
M = eqn.M(freeEdge,freeEdge); 
B = eqn.B(:,freeEdge);
f = eqn.f(freeEdge);
g = eqn.g;
u = zeros(length(freeEdge),1);

%% Isolate equations
pdeg = full(sum(spones(B),2));
isIsoElem = (pdeg == 1);
isoElem = find(isIsoElem);
[ii,jj,ss] = find(B(isIsoElem,:));
u(jj) = g(isoElem(ii))./ss;
temp = abs(g-B*u);
idx = temp<1e-14;
if any(idx)
    findelem(node,elem,idx,'index','FaceColor','c');  
end

%% Data structure
N = max(elem(:));
NT = size(elem,1);
NE = size(edge,1);
validEdge = 1:NE;
validEdge(freeEdge(jj)) = []; % remove isolate edges
e2v = sparse([validEdge,validEdge], double(edge(validEdge,:)), 1, NE, N);
validElem = find(~isIsoElem); % remove isolate element
t2v = sparse([validElem,validElem,validElem], elem(validElem,:), 1, NT, N);
oppe2v = sparse(double(elem2edge(:)), elem(:), 1, NE, N);
% In e2v, remove isolated edges dof since their values are determined. But
% they are remained in oppe2v as boundary dofs

%% Index map and indicator of elements
idxMap = zeros(NE,1);
idxMap(freeEdge) = 1:length(freeEdge);
isOpenElem = true(NT,1);
isSolvedElem = false(NT,1);
isSolvedElem(isoElem) = true; 
% In solved elements: Bu = g

%% Overlapping Schwarz smoother
scanNode = 1:N;
% scanNode = findboundary(elem);
% while ~isempty(scanNode)
for k = 1:length(scanNode)
    i = scanNode(k);
    starElem = find(t2v(:,i));
    if all(isSolvedElem(starElem)) % all elements are solved
        continue;
    end
    starEdge = idxMap(find(e2v(:,i))); %#ok<*FNDSB>
    starBdEdge = idxMap(find(oppe2v(:,i))); % edge opposite to a vertex
    domainStarEdgeIdx = (starEdge == 0); % remove fixed edge d.o.f.
    starEdge(domainStarEdgeIdx) = [];
    domainBdEdgeIdx = (starBdEdge == 0);
    starBdEdge(domainBdEdgeIdx) = [];
    % solve loccal problem
    locM = M(starEdge,[starEdge; starBdEdge]);
    locB = B(starElem,[starEdge; starBdEdge]);
    % update right hand side
    locf = f(starEdge) - locM*u([starEdge; starBdEdge]);
    locg = g(starElem) - locB*u([starEdge; starBdEdge]);
    locF = [locf; locg];
	Nf = length(locf);
    Ng = length(locg);
    locM = locM(:,1:Nf);
    locB = locB(:,1:Nf);
    locL = [locM locB'; locB  zeros(Ng)];
    % local problem is unique up to a constant
    % leave an open element non-solved
    % choice of open element: it is open and contains a max vertex number
    locOpenElem = find(isOpenElem(starElem));
    if ~isempty(locOpenElem)
        [maxN,I] = max(max(elem(starElem(locOpenElem),:),[],2));
        locOpenElem = locOpenElem(I);
    else
        [maxN,I] = max(max(elem(starElem,:),[],2));
        locOpenElem = I;
%         locOpenElem = Ng; % chose the last one if no elem open
    end
    activeIdx = true(Nf+Ng,1);
    activeIdx(Nf+locOpenElem) = false;
    locep = locL(activeIdx,activeIdx)\locF(activeIdx);
    u(starEdge) = u(starEdge) + locep(1:Nf);
    isSolvedElem(starElem) = true;
    locResidual = abs(locg(locOpenElem) - locB(locOpenElem,:)*locep(1:Nf));
    if locResidual > 1e-14 % check the residual of the open element
        isSolvedElem(starElem(locOpenElem)) = false;
    end
    isOpenElem(starElem) = false;
    %% debug: plot open elem and check residual of each element
    findelem(node,elem,starElem(locOpenElem),'noindex','FaceColor','r'); 
    findelem(node,elem,starElem(find(activeIdx(Nf+1:end))),'noindex','FaceColor','c'); 
%     temp = abs(g-B*u);
%     idx = temp<1e-14;
%     if any(idx)
%         findelem(node,elem,idx,'noindex','FaceColor','c');  
%     end
    findnode(node,i);
    pause(1)
end
% nextScanNode = false(N,1);
% nextScanNode(elem(~isOpenElem,:)) = true;
% nextScanNode(scanNode) = false;
% scanNode = find(nextScanNode);
% unSolvedElem = find(~isSolvedElem);
% scanNode = elem(unSolvedElem,:);
% scanNode = scanNode(:);
% end