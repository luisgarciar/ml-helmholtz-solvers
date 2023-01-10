%% P3 Cubic Element in 2D
%
% We explain degree of freedoms for the cubic element on triangles There
% are three types of dofs: vertex type, edge type and element type Given a
% mesh, the required data structure can be constructured by
%
%   [elem2dof,elem2edge,edge,bdDof,freeDof] = dofP3(elem)

%%
help dofP3

%% Local indexing of DOFs
node = [0 0; 1 0; 0.5 0.5*sqrt(3)];
elem = [1 2 3];
lambda = [1 0 0; 0 1 0; 0 0 1; ...  % 1,2,3 three vertices
          0 1/3 2/3; 0 2/3 1/3; ... % 4, 5 first edge
          2/3 0 1/3; 1/3 0 2/3; ... % 6, 7 second edge 
          1/3 2/3 0; 2/3 1/3 0; ... % 8, 9 third edge
          1/3 1/3 1/3];             % 10   center of element
dofNode = lambda*node;
figure;
set(gcf,'Units','normal'); 
set(gcf,'Position',[0,0,0.3,0.3]);
showmesh(node,elem);
hold on;
findnode(dofNode);

%% A local basis of P3
%
% The 10 Lagrange-type basis functions are denoted by $\phi_i, i=1:10$, i.e
% $\phi_i(x_j)=\delta _{ij},i,j=1:10$. In barycentric coordinates, they are
%
% $$ \phi_1 = 1/2(3\lambda_1-1) (3\lambda_1-2)\lambda_1,\quad \nabla \phi_1 = (27/2 \lambda_1 \lambda_1-9 \lambda_1+1) \nabla \lambda_1,$$
%
% $$ \phi_2 = 1/2(3\lambda_2-1) (3\lambda_2-2)\lambda_2,\quad  \nabla \phi_2 = (27/2 \lambda_2 \lambda_2-9 \lambda_2+1) \nabla \lambda_2,$$ 
%
% $$ \phi_3 = 1/2(3\lambda_3-1) (3\lambda_3-2)\lambda_3,\quad  \nabla \phi_3 = (27/2 \lambda_3 \lambda_3-9 \lambda_3+1) \nabla \lambda_3,$$ 
%
% $$ \phi_4 = 9/2\lambda_3\lambda_2 (3\lambda_2-1),\quad  \nabla\phi_4 = 9/2 ((3 \lambda_2 \lambda_2-\lambda_2) \nabla \lambda_3+ \lambda_3 (6 \lambda_2-1) \nabla \lambda_2)$$ 
%
% $$ \phi_5 = 9/2\lambda_3\lambda_2 (3\lambda_3-1),\quad  \nabla\phi_5 = 9/2 ((3 \lambda_3 \lambda_3-\lambda_3) \nabla \lambda_2+ \lambda_2 (6 \lambda_3-1) \nabla \lambda_3),$$ 
%
% $$ \phi_6 = 9/2\lambda_1\lambda_3 (3\lambda_3-1),\quad  \nabla\phi_6 = 9/2 ((3 \lambda_3 \lambda_3-\lambda_3) \nabla \lambda_1+ \lambda_1 (6 \lambda_3-1) \nabla \lambda_3),$$
% 
% $$ \phi_7 = 9/2\lambda_1\lambda_3 (3\lambda_1-1),\quad  \nabla\phi_7 = 9/2 ((3 \lambda_1 \lambda_1-\lambda_1) \nabla \lambda_3+ \lambda_3 (6 \lambda_1-1) \nabla \lambda_1),$$ 
%
% $$ \phi_8 = 9/2\lambda_1\lambda_2 (3\lambda_1-1),\quad  \nabla\phi_8 = 9/2 ((3 \lambda_1 \lambda_1-\lambda_1) \nabla \lambda_2+ \lambda_2 (6 \lambda_1-1) \nabla \lambda_1),$$ 
%
% $$ \phi_9 = 9/2\lambda_1\lambda_2 (3\lambda_2-1),\quad  \nabla\phi_9 = 9/2 ((3 \lambda_2 \lambda_2-\lambda_2) \nabla \lambda_1+ \lambda_1 (6 \lambda_2-1) \nabla \lambda_2),$$ 
%
% $$ \phi_{10} = 27\lambda_1\lambda_2\lambda_3, \quad \nabla \phi_{10} =  27 (\lambda_1 \lambda_2 \nabla \lambda_3+\lambda_1 \lambda_3 \nabla \lambda_2+ \lambda_3 \lambda_2 \nabla \lambda_1).$$
%
% When transfer to the reference triangle formed by $(0,0),(1,0),(0,1)$,
% the local bases in x-y coordinate can be obtained by substituting 
% 
% $$\lambda _1 = x, \quad \lambda _2 = y, \quad \lambda _3 = 1-x-y.$$ 

%% Global indexing of DOFs
%
%  The global indices of the dof is organized  according to the order of
%  nodes, edges and elements. To be consistent, the dof on an edge depends
%  on the orientation of edge only. 
