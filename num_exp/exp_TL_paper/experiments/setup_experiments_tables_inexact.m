kk = [10 20 40 80 100];

%factoreps = 0.5;
%p = gmres2D_csl_vs_dcsl_kvarying_coarse_inexact(kk,factoreps);

factoreps = 2;
p = gmres2D_csl_vs_dcsl_kvarying_coarse_inexact(kk,factoreps);

factoreps = 5;
p = gmres2D_csl_vs_dcsl_kvarying_coarse_inexact(kk,factoreps);
