%Test multigrid Sommerfeld BC's

clc;
bc  = 'som'
dim = 2;
npf = 5;

R = fwrestriction_som(npf,dim,bc);

R=16*full(R);
          
R(1,:)
