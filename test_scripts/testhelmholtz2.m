%script for testing the function helmholtz2.m

k  = 10;
b1 = 1;
b2 = 0.5;
np = ceil( 10 * k / pi);
bc = 'som'

A = helmholtz2(k,np,np,bc);
S = shift_laplace2(k,b1,b2,np,np,bc);

egv = eig(full(S\A));
%egv = eig(full(A));

plot(real(egv),imag(egv),'r*');
axis equal