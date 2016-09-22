%Test multigrid Sommerfeld BC's
clc;
bc  = 'som';
dim = 2;
k   = 50;
npc = 50; npf = 2*(npc)+1; eps=0.5*k^2;

A = helmholtz2(k,eps,npf,npf,bc);
tic
restr  =  fwrestriction_som(npf,dim,bc); %size(restr)
interp =  4*restr';
time1  = toc;

%Two-grid cycle
npre = 2; npos = 2; numcycles = 10;
b       = zeros(length(A),1); %b(ceil(length(A)/2),1)=1;
x_init  = randn(length(A),1);

L=sparse(tril(A));
U=sparse(triu(A)); 
D=sparse(diag(diag(A)));

%gauss-seidel smoothing by default
tic
[x_sol,relres] = twogrid_som(A,L,U,D,restr,interp,b,x_init,npre,npos,numcycles);
time2=toc;
