lev = 5;                        %number of levels
npc = 1;                        %number of points in coarsest grid (1D) 
npf = (2^lev)*npc+(2^lev-1);    %number of points in finest grid (1D)
bc  = 'dir';
b1  = 1; b2 = 0.5;
k   = 10;
dim = 1;

npc = 1; npf = 2*(npc)+1;
interp   = fwrestriction(npf,dim); size(interp)
restr    = lininterpol(npc, dim);   size(restr)

%[SLgrid_matrices,restrict,interp] = mg_setup(M,lev,bc,dim);
