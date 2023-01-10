%% P1 Linear Element
%
% For the linear element on a simplex, the local basis functions are
% barycentric coordinate of vertices. The local to global pointer is
% |elem|. This is the default element for elliptic equations.

%% A local basis of P1
% 
% For $i = 1, 2,..., d+1$, the local basis of linear element space is
%
% $$\phi_i = \lambda_i, \nabla \phi_i = \nabla \lambda_i = - \frac{|e_i|}{d!|T|}\mathbf n_i,$$
%
% where $e_i$ is the edge opposite to the i-th vertex and $n_i$ is the unit
% outwards normal direction.
%
% See <http://www.math.uci.edu/~chenlong/226/Ch2FEM.pdf Finite Element
% Methods> Section 2.1 for geometric explanation of the barycentric
% coordinate.

%% Global indexing of DOFs
node = [0,0; 1,0; 1,1; 0,1];
elem = [2,3,1; 4,1,3];      
[node,elem] = uniformbisect(node,elem);
figure;
showmesh(node,elem);
findnode(node);
findelem(node,elem);
display(elem);