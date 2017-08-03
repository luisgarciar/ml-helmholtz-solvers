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
bdFlag = setboundary(node,elem,'Robin'); 
