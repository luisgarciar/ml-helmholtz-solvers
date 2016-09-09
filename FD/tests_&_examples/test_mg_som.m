%Test multigrid Sommerfeld BC's

clc;
bc  = 'som'
dim = 2
npi = 3
R = full(fwrestriction_som(npi,dim,bc))
size(R)

