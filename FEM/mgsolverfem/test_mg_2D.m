%Multigrid SL FEM test

k = 10;
poweps    = 1;
factoreps = 1;
bc = 'som';

pollution = 'no';
npf = ceil(k^(3/2));
    
[node,elem] = squaremesh([0 1 0 1],0.5);
   for j = 1:9
     [node,elem] = uniformrefine(node,elem);
   end
   
[bdNode,bdEdge,isBdNode] = findboundary(elem);
bdFlag  = setboundary(node,elem,'ABC');
pdeSL   = helmholtz2Dconstantwndata(k,factoreps,poweps);

[eqn2,~] = helmholtz2Dfem(node,elem,pdeSL,bdFlag,bdEdge);
    
%shifted Laplace matrix
Aeps = eqn2.A;
n = size(Aeps,1);
b = sparse(n,1);
option.solver = 'NO';

[x,info,Ai,Bi,BTi,Res,Pro,isFreeDof] = mg(Aeps,b,elem,option);

[HB,NL,level] = HBstructure(elem,50);

