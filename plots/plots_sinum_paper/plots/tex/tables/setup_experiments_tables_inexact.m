kk = [10 20 30 40];

factoreps = 0.5;
p = gmres2D_csl_vs_dcsl_kvarying_coarse_inexact(kk,factoreps);

%factoreps = 1;
%p = gmres2D_csl_vs_dcsl_kvarying_coarse_inexact(kk,factoreps);

%factoreps = 2;
%p = gmres2D_csl_vs_dcsl_kvarying_coarse_inexact(kk,factoreps);
