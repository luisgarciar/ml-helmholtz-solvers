%Test for function slapfov.m
%clear all
%parameters of Helmholtz problem and shifted Laplacian
k         = 100;
poweps    = 2;
factoreps = 1;
bc = 'som';
np = ceil(k^(3/2));

fvpts = 50;%number of points for fov plot

if (mod(np+1,2)==1) 
    np = np+1; 
end
npc = (np-1)/2;

A2    = helmholtzfem(k,np,0,bc);   %Helmholtz matrix
Aeps2 = helmholtzfem(k,np,eps,bc); %Shifted Laplace matrix

[fovAhat2,minfov] = slapfov(A2,Aeps2,fvpts); %field of values (complex valued vector)  

reFOV2 = real(fovAhat2); imFOV2 = imag(fovAhat2);
cvh   = convhull(reFOV,imFOV);
 
plot(reFOV2(cvh),imFOV2(cvh),'LineWidth',2)
hold on
%plot(0,0,'k+','Markersize',6,'LineWidth',2)
 %axis equal