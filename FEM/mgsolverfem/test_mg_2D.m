 k  = 10;
 poweps    = 1;
 factoreps = 1;
 bc = 'som';
 
% %coarse mesh
 npcc = 10; 
 h = 1/(npcc+1);
 numlev = 5;
 
pdeSL = helmholtz2Dconstantwndata(k,factoreps,poweps);

%option.twolevel = true;
option.twolevel = false;
[mg_mat,mg_split,restr,interp] = mg_setupfem_2D(npcc,numlev,pdeSL,option);

% [node,elem] = squaremesh([0 1 0 1],h);
% 
% %refining the mesh numlev times
% for j = 1:numlev-1
%     [node,elem] = uniformrefine(node,elem);
% end
% 
% %Find boundary nodes
% [~,bdEdge,~] = findboundary(elem);
%     
% %Sets Sommerfeld boundary conditions on all boundary edges
% bdFlag = setboundary(node,elem,'ABC');
% 
% %pde data
% pde     = helmholtz2Dconstantwndata(k,factoreps,poweps);
% [eqn,~] = helmholtz2Dfem(node,elem,pde,bdFlag,bdEdge);
% 
% A = eqn.A;
% n = size(A,1);
% b = sparse(n,1);
% option.solver = 'NO';
% 
% [~,~,Ai,~,~,Res,Pro,~] = mg(A,b,elem,option);
% assert(length(Ai)==numlev, 'error: incorrect number of levels');
% 
% mg_mat = flip(Ai);
% restr  = flip(Res);
% restr  = restr(1:numlev-1);
% 
% prol  = flip(Pro); 
% prol  = prol(2:numlev);
% 
% mg_split = cell(numlev,1);    
%
% for i=1:numlev     
%     %matrix splitting of mg_mat{i}
%     mg_split{i}.U = sparse(triu(mg_mat{i},1));  
%     mg_split{i}.L = sparse(tril(mg_mat{i},-1));
%     mg_split{i}.D = spdiags(diag(mg_mat{i}),0,length(mg_mat{i}),length(mg_mat{i}));
% 
% end

