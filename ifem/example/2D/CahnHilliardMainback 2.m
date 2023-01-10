function CahnHilliardMain
%% 
%\frac{\partial c}{\partial t} - \nabla(\nabla(\frac{df}{dc}-\lambda\nabla^{2}c))= 0 \\Omega,
%(\nabla\left(\frac{d f}{d c} - \lambda \nabla^{2}c\right)\right)= 0 \partial\Omega
%\lambda \nabla c \cdot n &=& 0 \partial\Omega.

%  http://fenicsproject.org/documentation/dolfin/dev/python/demo/pde/cahn-hilliard/python/documentation.html
close all
%% Parameters
%  maxN = 2e3;     theta = 0.5;    maxIt = 1; 
%  N = zeros(maxIt,1);     errL2 = zeros(maxIt,1);     errH1 = zeros(maxIt,1);

%%  Generate an initial mesh
% node = [0,0; 0,1; 1,0;1,1];            % nodes
% elem = [1,3,2; 4,2,3];                 % elements
% elem = label(node,elem);               % label the mesh
h   = 1/100;
[node,elem] = squaremesh([0,1,0,1],h);

%bdEdge = setboundary(node,elem,'Neumann');       % Dirichlet boundary condition    
% showmesh(node,elem);                            % plot mesh                
% findelem(node,elem,'all','index','color','g');  % plot element indices
% findnode(node,'all','index','color','r');       % plot node indices
%%  Get a fine mesh by uniform bisection
% for k = 1:0
%      [node,elem,bdEdge] = uniformbisect(node,elem,bdEdge);  
% %    [node,elem] = bisect(node,elem,'all');
% end

% showmesh(node,elem);                            % plot mesh                
% findelem(node,elem,'all','index','color','g');  % plot element indices
% findnode(node,'all','index','color','r');       % plot node indices
% NT  = size(elem,1);
N   = size(node,1);
%% Parametere
global efsilonsquare 
efsilonsquare = 1.0e-2;
pde = CahnHilliarddata(h);
tao  = 5.0e-6;  %time step
times = 80;     
%% initial condition;
u0 = pde.initial_u(node); % initial vavlue u
w0 = pde.initial_w(node); % initial vavlue w
%% Solve
<<<<<<< mine
[tempvar,tempvar,u0] = CahnHilliardP1(node,elem,pde,w0,u0,tao,times);
=======
[~,~,u0] = CahnHilliardP1(node,elem,pde,w0,u0,tao,times);
%[~,u0] = CahnHilliardP2(node,elem,pde,w0,u0,tao,times);
>>>>>>> theirs
%% Plot
x   = node(:,1);
y   = node(:,2);
z   = full(u0);
x0   = reshape(x,sqrt(N),sqrt(N));
y0   = reshape(y,sqrt(N),sqrt(N));
z0   = reshape(z,sqrt(N),sqrt(N));
surf(x0,y0,z0);
% shading interp
% [X,Y,Z]=griddata(x,y,z,linspace(min(x),max(x),100)',linspace(min(y),max(y),100),'v4');%插值
% figure,contourf(X,Y,Z); %等高线图
%colorbar



