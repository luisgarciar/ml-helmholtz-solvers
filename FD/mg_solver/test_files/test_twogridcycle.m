bc  = 'dir';
b1  = 1; b2 = 0.5;
k   = 10;
dim = 2;

npc = 50; npf = 2*(npc)+1;

switch dim
    case  1
        A = helmholtz(k,npf,npf,bc);
    case 2
        A = helmholtz2(k,npf,npf,bc);
end

restr  =  fwrestriction(npf,dim); %size(restr)
interp =  lininterpol(npc, dim); 


%%
npre = 2; npost = 2; smo = 'wjac'; w = 2/3; numcycles = 11;

b       = zeros(length(A),1); %b(ceil(length(A)/2),1)=1;
x_init  = randn(length(A),1);

profile clear
profile on
[x_sol,relres] = twogrid(A,restr,interp,b,x_init,dim,npre,npost,w,smo,numcycles);
profile off