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

kk  =  [20 40 80 160 320];

%CSL
opts.prec = 'csl';
[fovAAeps] = computefovAhat(kk,opts);
plotfovAhat(fovAAeps,kk,opts);
close all

%Deflated CSL
opts.prec = 'adef';
[fovAB] = computefovAhat(kk,opts);
plotfovAhat(fovAB,kk,opts);
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

%kmult  =  [20 40 80 160 320];

%CSL
opts.prec = 'csl';
[fovAAeps] = computefovAhat(kk,opts);
plotfovAhat(fovAAeps,kk,opts);
close all

%Deflated CSL
opts.prec = 'adef';
[fovAB] = computefovAhat(kk,opts);
plotfovAhat(fovAB,kk,opts);
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

kk  =  [20 40 80 160 320];

%CSL
opts.prec = 'csl';
[fovAAeps] = computefovAhat(kk,opts);
plotfovAhat(fovAAeps,kk,opts);
close all

%Deflated CSL
opts.prec = 'adef';
[fovAB] = computefovAhat(kk,opts);
plotfovAhat(fovAB,kk,opts);
close all

%% Plots with  N = ceil(k^1.5)
%1D Case
opts.dim = 1;
opts.poweps    = 2;
opts.factoreps = 5;
opts.bc = 'som';
opts.fvpts = 64;
opts.disc = 'pf';

%kk  =  [20 40 80 160 320];

%CSL
opts.prec = 'csl';
[fovAAeps] = computefovAhat(kk,opts);
plotfovAhat(fovAAeps,kk,opts);
close all

%Deflated CSL
opts.prec = 'adef';
[fovAB] = computefovAhat(kk,opts);
plotfovAhat(fovAB,kk,opts);
close all
