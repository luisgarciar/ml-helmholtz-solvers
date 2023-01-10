k  = 20;
poweps    = 2;
factoreps = 1;
bc = 'som';
par = 0.5;
npcc = 3;

[npf,numlev] = fem_npc_to_npf(npcc,k,par);  %number of points in finest grid (1D)
h= 1/(npf+1);

pdeSL = helmholtz2Dconstantwndata(k,factoreps,poweps);
option.twolevel = false;

[mg_mat,mg_split,restr,interp] = mg_setupfem_2D(npcc,numlev,pdeSL,option);
[node,elem] = squaremesh([0 1 0 1],h);

Aeps = mg_mat{1};
N = length(Aeps);
b=zeros(N,1);
b(ceil(N/2),1)=1;

[L,U] = lu(Aeps);
u_gal = U\(L\b);

npre = 1; npos = 1; w  = 0.7; smo = 'wjac';
u0      = zeros(length(Aeps),1);
Aepsinv = @(x,num) Wcycle(mg_mat,mg_split,restr,interp,u0,x,npre,npos,w,smo,num);
u_mg = Aepsinv(b,10);

%%
figure(1)
showsolution(node,elem,real(full(u_gal)));

figure(2)
showsolution(node,elem,real(full(u_mg)));

norm(u_mg-u_gal,Inf)




