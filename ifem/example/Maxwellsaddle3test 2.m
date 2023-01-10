close all; clear all
pde = Maxwellsaddledata;
% cube
% [node,elem] = cubemesh([0,1,0,1,0,1],0.25);
% % bdFlag = setboundary3(node,elem,'Neumann'); % Pure Neumann boundary condition doesn't work.
% bdFlag = setboundary3(node,elem,'Dirichlet');
% Lshape
[node,elem] = cubemesh([-1,1,-1,1,-1,1],0.5);
% [node,elem] = delmesh(node,elem,'x>0 & y<0');
[node,elem] = delmesh(node,elem,'x<0 & y<0 & z>0');
% showboundary3(node,elem);
bdFlag = setboundary3(node,elem,'Dirichlet');

err = zeros(4,1); 
N = zeros(4,1);
% option.solver = 'diag';
option.solver = 'tri';
% option.solver = 'direct';
option.printlevel = 2;
option.mg.Vit = 1;
option.mg.smoothingstep = 3;    % Smoothing step.
option.mg.smoothingratio = 1.5; % ratio of variable smoothing
% option.solver = 'direct';
for i = 1:4
    u = Maxwellsaddle(node,elem,pde,bdFlag,option);
    err(i) = getHcurlerror3NE(node,elem,pde.curlu,real(u));
    N(i) = size(u,1);
    [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
end
showrate(N,err,2);