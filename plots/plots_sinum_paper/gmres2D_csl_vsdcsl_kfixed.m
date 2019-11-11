%Comparison of Shifted Laplacian and two-level deflation for fixed k
%Homogeneous problem



dim      = 2;
poweps   = 2;
factoreps = 1;

kk = 10*pi;
bc = 'som';

reflevs = 1;   %
restart   = [];
tol       = 1e-8;
maxit     = 100;


npc = ceil(k^(3/2)/8);
npf = (2*npc+1);


 H = 1/(npc+1);
 [node,elem] = squaremesh([0,1,0,1],H);  %coarse mesh
 
 for i=1:reflevs
     [node,elem] = uniformrefine(node,elem); %fine mesh (uniformly refined)
 end
 
%set boundary conditions 
[bdNode,bdEdge,isBdNode] = findboundary(elem);

