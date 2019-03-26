%Plots of FOVs for SINUM paper

%% Plots with eps = k^2

%% Plots with N = ceil(k^2)
%1D Case
opts.dim = 1;
opts.poweps    = 2;
opts.factoreps = 1;
opts.bc = 'som';
opts.fvpts = 64;
opts.disc = 'q'; 

kmult  =  [5 20 50 80 100];

%CSL
opts.prec = 'csl';
[fovAAeps] = computefovAhat(kmult,opts);
plotfovAhat(fovAAeps,kmult,opts);
close all

%Deflated CSL
opts.prec = 'adef';
[fovAB] = computefovAhat(kmult,opts);
plotfovAhat(fovAB,kmult,opts);
close all

%% Plots with  N = ceil(k^1.5)
%1D Case
opts.dim = 1;
opts.poweps    = 2;
opts.factoreps = 1;
opts.bc = 'som';
opts.fvpts = 64;
opts.disc = 'pf';
close all

kmult  =  [5 20 50 80 100];

%CSL
opts.prec = 'csl';
[fovAAeps] = computefovAhat(kmult,opts);
plotfovAhat(fovAAeps,kmult,opts);
close all

%Deflated CSL
opts.prec = 'adef';
[fovAB] = computefovAhat(kmult,opts);
plotfovAhat(fovAB,kmult,opts);
close all




%% Plots with eps = 5k^2

%% Plots with N = ceil(k^2)
%1D Case
opts.dim = 1;
opts.poweps    = 2;
opts.factoreps = 5;
opts.bc = 'som';
opts.fvpts = 64;
opts.disc = 'q';

kmult  =  [5 20 50 80 100];

%CSL
opts.prec = 'csl';
[fovAAeps] = computefovAhat(kmult,opts);
plotfovAhat(fovAAeps,kmult,opts);
close all

%Deflated CSL
opts.prec = 'adef';
[fovAB] = computefovAhat(kmult,opts);
plotfovAhat(fovAB,kmult,opts);
close all



%% Plots with  N = ceil(k^1.5)
%1D Case
opts.dim = 1;
opts.poweps    = 2;
opts.factoreps = 5;
opts.bc = 'som';
opts.fvpts = 64;
opts.disc = 'pf';

kmult  =  [5 20 50 80 100];

%CSL
opts.prec = 'csl';
[fovAAeps] = computefovAhat(kmult,opts);
plotfovAhat(fovAAeps,kmult,opts);
close all

%Deflated CSL
opts.prec = 'adef';
[fovAB] = computefovAhat(kmult,opts);
plotfovAhat(fovAB,kmult,opts);
close all


%% Plots with eps = 10k^2

%% Plots with N = ceil(k^2)
%1D Case
opts.dim = 1;
opts.poweps    = 2;
opts.factoreps = 10;
opts.bc = 'som';
opts.fvpts = 64;
opts.disc = 'q';

kmult  =  [5 20 50 80 100];

%CSL
opts.prec = 'csl';
[fovAAeps] = computefovAhat(kmult,opts);
plotfovAhat(fovAAeps,kmult,opts);
close all

%Deflated CSL
opts.prec = 'adef';
[fovAB] = computefovAhat(kmult,opts);
plotfovAhat(fovAB,kmult,opts);
close all


%% Plots with  N = ceil(k^1.5)
%1D Case
opts.dim = 1;
opts.poweps    = 2;
opts.factoreps = 10;
opts.bc = 'som';
opts.fvpts = 64;
opts.disc = 'pf';

kmult  =  [5 20 50 80 100];

%CSL
opts.prec = 'csl';
[fovAAeps] = computefovAhat(kmult,opts);
plotfovAhat(fovAAeps,kmult,opts);
close all

%Deflated CSL
opts.prec = 'adef';
[fovAB] = computefovAhat(kmult,opts);
plotfovAhat(fovAB,kmult,opts);
close all


