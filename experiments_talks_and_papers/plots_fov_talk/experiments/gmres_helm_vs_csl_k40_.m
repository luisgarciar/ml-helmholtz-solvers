
%% Solution of 2D Helmholtz equation for k=20 and point source

%% Construction of the matrices
clear all
close all
save_flag = 1;  % save_flag=1: save plots and table, =0 do not save.

% Setup parameters
% Setup list of wavenumbers and shifts


k  = 40; %Wavenumber
bc = 'som'; % boundary conditions
npf = ceil(k^(3/2)); %number of gridpoints (no pollution)

if (mod(npf+1,2)==1)  %set an even number of interior points in 1D
    npf = npf+1;
end

%construction of the mesh
h = 1/(npf+1);
[node,elem] = squaremesh([0,1,0,1],h);  %coarse mesh
[bdNode,bdEdge,isBdNode] = findboundary(elem);

bdFlag = setboundary(node,elem,'ABC');
pdehelm = helmholtz2Dconstantwndata(k,0,1);

option.tol = 1e-8;
[eqn1,~] = helmholtz2Dfem(node,elem,pdehelm,bdFlag,bdEdge);

A = eqn1.A;
f = zeros(length(A),1); f(ceil(length(f)/2))=1;
u = real(A\f);

showsolution(node,elem,u);
box off
set(findobj(gcf, 'type','axes'), 'Visible','off')
%matlab2tikz('helmholtzk20.tex','standalone', true);





