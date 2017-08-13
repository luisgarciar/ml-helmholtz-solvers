
k   = 5;
dim = 2;
pollution = 'no';
npf = ceil(k^(3/2));
        
    if (mod(npf+1,2)==0)  %set an odd number of interior points in 1D
        npf = npf+1;
    end
    npc = (npf-1)/2;
    

scale = 2.^(2:2:8);    


    for jj = 1:length(scale)

    poweps    = 2;
    factoreps = 1;
    bc = 'som';

    %Construct square mesh of meshsize h
    h = 1/npf;
    [node,elem] = squaremesh([0,1,0,1],h);

    %Visualize mesh, nodes, elements
    showmesh(node,elem)
    hold on;
    findnode(node,1:length(node));% find node indices

    %Find boundary nodes
    [bdNode,bdEdge,isBdNode] = findboundary(elem); 

    %Sets Sommerfeld boundary conditions on all boundary edges
    bdFlag = setboundary(node,elem,'ABC'); 

    %the structure pde contains data for a simple test problem
    pde = helmholtz2Dtrigdata(k);

    [eqn,info] = helmholtz2Dfem(node,elem,pde,bdFlag);

    %Matrix and right hand side
    A=eqn.A; b=eqn.b;
    %b = zeros(length(A),1); b(ceil(length(A)/2),1)=1;
    u = A\b;
    u_exact = pde.exactu(node);

    showsolution(node,elem,u);
    end
