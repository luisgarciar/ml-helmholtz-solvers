%Minimum eigenvalue of Helmholtz matrices

kk = 20:400;
%kk =10;

N= length(kk);
mineigv = zeros(N,1);
taneigv = zeros(N,1);

for i=1:N
    k=kk(i);
    n = ceil(k/0.625);
    %n = ceil(k^(3/2));
    h = 1/n;
    j = 1:(n-1);
    eigvA = 4*sin(j*pi*h/2).^2-k^2*h^2;
    mineigv(i) = min(abs(eigvA));
    taneigv(i) = k/min(abs(eigvA));   
end

semilogy(kk,mineigv);
hold on
semilogy(kk,1./kk,'r');

