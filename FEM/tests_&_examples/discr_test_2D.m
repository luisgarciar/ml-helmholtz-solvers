%Test of discretization of 2D Helmholtz problem with P1 finite elements
%
% We use the test problem
%
%  -div(grad u)-k^2 u   = 0   in Omega= (0,1)x(0,1)
%  grad(u) dot n - i*ku = g on boundary(Omega) 
%
%  where the function g is chosen such that
%  the exact solution is a planewave, i.e.,
%  u = e^{i kk \cdot (x,y)}
%  where kk = k (cos(t), sin(t)) 
%
% Reference: 
% Finite Element Analysis of Acoustic Scattering,
% Ihlenburg (1997) - 
% Section 1.1, pp. 5 and section 4.2.2, p. 108

%Set wavenumber and number of points
k   = 2*pi;
npf = ceil(k^(3/2));

%np = npf;
np = npf*2.^(2:2:6);
plot_sol = 'no';

relerror = zeros(length(np),1);

for i=1:length(np)
    npf = np(i);
    if (mod(npf+1,2)==0)  %set an odd number of interior points in 1D
        npf = npf+1;
    end
    npc = (npf-1)/2;
    
    %Construct square mesh of meshsize h
    h = 1/npf;
    [node,elem] = squaremesh([0,1,0,1],h);
    
    %Find boundary nodes
    [bdNode,bdEdge,isBdNode] = findboundary(elem);
    
    %Set Sommerfeld boundary conditions on all boundary edges
    bdFlag = setboundary(node,elem,'ABC');
    
    %the structure pde contains data for the planewave test problem
    t   = pi/2; %propagation direction
    pde = helmholtz2Dplanewavedata(k,t);
    option.tol = 1e-12;
    
    %Constructing the matrix and right hand side
    [eqn,info] = helmholtz2Dfem(node,elem,pde,bdFlag,bdEdge);
   
    
    %Matrix and right hand side
    A = eqn.A; b = eqn.b;
    u       = A\b;
    
    %interpolant of the exact solution
    u_exact = pde.exactu(node);

    %comparison of the two solutions (exact vs. computed)
    
    if strcmp(plot_sol,'yes')
        figure(1)
        showsolution(node,elem,real(u_exact));
        figure(2)
        showsolution(node,elem,real(u));
       
    end
  
    %relative error in euclidean norm;
    relerror(i) = norm(u-u_exact)/norm(real(u_exact));
end


relerror