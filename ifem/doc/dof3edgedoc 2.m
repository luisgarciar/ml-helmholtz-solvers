%% Dof on Edges in Three Dimensions
%
% We describe general idea of the data structures generated in subroutine 
% dof3edge for three dimensional tetrahedron grids. 

help dof3edge

%% Local Labeling of DOFs
node = [1,0,0; 0,1,0; 0,0,0; 0,0,1];
elem = [1 2 3 4];
localEdge = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
edge = zeros(20,2);
edge([1 12 5 20 11 4],:) = localEdge;
elem2edge = [1 12 5 20 11 4];
figure(1); clf;
set(gcf,'Units','normal'); 
set(gcf,'Position',[0,0,0.6,0.4]);
subplot(1,2,1)
showmesh3(node,elem);
view(-14,12);
findnode3(node);
findedge(node,localEdge,'all');
subplot(1,2,2)
showmesh3(node,elem);
view(-14,12);
% findnode3(node);
findedge(node,edge,elem2edge');
%%
% The six dofs associated to edges in a tetrahedron is sorted in the
% ordering [1 2; 1 3; 1 4; 2 3; 2 4; 3 4]. Here [1 2 3 4] are local indices
% of vertices.
%
% For a triangulation consisting of several tetrahedrons, all local edges
% are collected and the duplication is eliminated to form edges of the
% triangulation. Therefore the global indices of edges is different with
% the local one. elem2edge(:,1:6) records the mapping from local to global
% indices.
display(elem2edge)

%% Orientation of Edges
node = [1,0,0; 0,1,0; 0,0,0; 0,0,1];
elem = [2 4 3 1];
localEdge = [elem(:,1) elem(:,2); elem(:,1) elem(:,3); elem(:,1) elem(:,4); ...
             elem(:,2) elem(:,3); elem(:,2) elem(:,4); elem(:,3) elem(:,4)];
edge = sort(localEdge,2);
elem2edge = [1 2 3 4 5 6];
elem2edgeSign = [1 1 -1 -1 -1 -1];     
figure(1); clf;
set(gcf,'Units','normal'); 
set(gcf,'Position',[0,0,0.6,0.4]);
subplot(1,2,1)
showmesh3(node,elem);
view(-14,12);
tempnode = node(elem(:),:);
findnode3(tempnode);
findedge(node,localEdge,'all','vec');
subplot(1,2,2)
showmesh3(node,elem);
view(-14,12);
findnode3(node);
findedge(node,edge,'all','vec');
%%
% The edge [i,j] is orientated in the direction such that i<j, where i,j
% are *global* indices. The edges formed by local indices may not be
% consistent with this orientation. Therefore elem2edgeSign is used to indicate
% the consistency.
%
% The nodal indices in the left figure is local while that in the right is
% the global one. The local direction and global direction of edges is also
% indicated by different edge vectors.
display(elem2edge)
display(elem2edgeSign)

%% An Example
[node,elem] = cubemesh([0,1,0,1,0,1],1);
[elem2edge,edge,elem2edgeSign] = dof3edge(elem);
figure(1); clf;
showmesh3(node,elem,[],'FaceAlpha',0.25);
view(-30,10);
findedge(node,edge,'all','vec');
findelem3(node,elem,[1 3]');
display(elem2edge);
display(elem2edgeSign);
