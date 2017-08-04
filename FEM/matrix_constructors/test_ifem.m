%Construct square mesh of meshsize h
h = 0.5;
[node,elem] = squaremesh([0,1,0,1],h);

%Visualize mesh, nodes, elements
showmesh(node,elem)
hold on;
findelem(node,elem); % plot element indices
findnode(node,1:length(node));% find node indices

%Find boundary nodes
[bdNode,bdEdge,isBdNode] = findboundary(elem); 

%Sets Robin boundary conditions on all boundary edges
bdFlag = setboundary(node,elem,'ABC'); 
bdFlag

ABC = [];
    isABC = (bdFlag(:) == 3);
    if any(isABC)
        allEdge = [elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])];
        ABC = allEdge(isABC,:);
    end

k=10;
ve = node(ABC(:,1),:) - node(ABC(:,2),:);
edgeLength = sqrt(sum(ve.^2,2)); 
ii = [ABC(:,1),ABC(:,1),ABC(:,2),ABC(:,2)];
jj = [ABC(:,1),ABC(:,2),ABC(:,1),ABC(:,2)];
temp = -sqrt(-1)*k*edgeLength;

%these 1/3,1/6 result from computing the boundary integrals

%int_{\Gamma} phi_i^2 ds     = 1/3*edgeLength. 
%int_{\Gamma} phi_i phi_j ds = 1/6*edgeLength.

ss = [1/3*temp, 1/6*temp, 1/6*temp, 1/3*temp]; 


      


