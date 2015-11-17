clear all


k     = 15
ppw   = 20
n     = 2^3;    %number of interior gridpoints
sigma = k^2;      %wavenumber

x     = [0:1/n:1]'; %1-D grid
x_res = x(2:length(x)-1);

f            = zeros(size(x));  %rhs function
f((n/2)+1,1) = 10; 
f            = f(2:length(f)-1);   

[A, sol] = helmholtz_1D(f,sigma,1);
[SLP,~ ] = helmholtz_1D(f,sigma*(1-0.5*1j),1) 


Id = eye(length(A));
m  = (length(A)+1)/2-1;

invSLP = inv(SLP)

Ih =zeros(m,length(A));

for i=1:length(A)
    Ih(:,i)=fwrestriction(Id(:,i));
end

save('file1.mat','A','invSLP','Ih')



