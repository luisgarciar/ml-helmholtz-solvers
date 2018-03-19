%Residual Comparison GMRES Helmholtz vs CSL  - 2D Problem

clear all;
save_flag = 1;  % save_flag=1: save plots and table, =0 do not save.

%Setup parameters
%Setup list of wavenumbers and shifts
%wavenum   = 20:20:100 ;
%poweps    = [1 1.5 2];
wavenum =  40;
poweps  = 1;



eps = k^poweps;
