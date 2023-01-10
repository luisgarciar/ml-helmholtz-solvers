%% FOV plot of Grcar Matrix


A = gallery('grcar',100);

fvA = fv(A,1,64);
eigvA = eig(full(A));


plot(real(fvA),imag(fvA),'b','Linewidth',3);
hold on
plot(real(eigvA),imag(eigvA),'k+');
axis equal
set(gca,'Xtick',[-2 0 2 4],'FontSize',35);
set(gca,'Ytick',[-2 0 2],'FontSize',35);

matlab2tikz('fovgrcar.tex','standalone', true);
