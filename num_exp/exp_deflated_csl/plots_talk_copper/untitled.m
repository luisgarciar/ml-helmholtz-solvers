% Plots for presentation
wavenum = 20; %% run this when testing changes in the code
pmin    = 12;       %min points per wavelenght
beta    = 0.5;      %complex shift
bc      = 'som' ; 

%  Sommerfeld Problem, spectrum k=10
A = helmholtz2(k,npx,npy,bc);