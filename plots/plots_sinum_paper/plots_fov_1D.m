%Plots of FOVs for SINUM paper

%1D Case
opts.dim = 1;
opts.poweps    = 2;
opts.factoreps = 1;
opts.bc = 'som';
opts.fvpts = 64;
opts.disc = 'pf';

kmult  =  [10 20 50 80 100];

%CSL
opts.prec = 'csl';
[fovAAeps] = computefovAhat(kmult,opts);


%Deflated CSL
opts.prec = 'adef';
[fovAB] = computefovAhat(kmult,opts);


