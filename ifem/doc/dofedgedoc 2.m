%% Dof on Edges in Two Dimensions
% 
% We describe general idea of the data structures generated in subroutine 
% dofedge for two dimensional triangular grids. 

help dofedge

%% Example
% Generate edge and elem2dof
[node,elem] = squaremesh([0,1,0,1],1/2);
T = auxstructure(elem);
elem2edge = T.elem2edge;
edge = T.edge;
edge2elem = T.edge2elem;

showmesh(node,elem);
findnode(node,'all');
findelem(node,elem,'all');
findedge(node,edge,'all','vec');

display(elem);
display(elem2edge);
display(edge);

% Consistency of oritentation of edges
NT = size(elem,1); NE = size(edge,1);
elem2edgeSign = ones(NT,3);
totalEdge = uint32([elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])]);
idx = (totalEdge(:,1)>totalEdge(:,2));
elem2edgeSign(idx) = -1;

display(elem2edgeSign);