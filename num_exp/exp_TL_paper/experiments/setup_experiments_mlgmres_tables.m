
factoreps = 1;
kk = [20 40 60 80 120 160];
p = mlgmres2D_csl_vs_dcsl_kvarying(kk,factoreps);

%factoreps = 1;
%p = gmres2D_csl_vs_dcsl_kvarying_coarse(kk,factoreps);

factoreps = 2;
p = mlgmres2D_csl_vs_dcsl_kvarying(kk,factoreps);
