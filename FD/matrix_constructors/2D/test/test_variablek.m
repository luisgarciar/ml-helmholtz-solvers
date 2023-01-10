%% Script for testing the function helmholtz2var.m

kmax  = 20; kmin = 30;    %For problem with random wavenumbers
kref  = 40;               %For wedge problem
np    = ceil(30*kref/pi); % number of grid points
npx   = np;
npy   = np;
hx    = 1/(npx+1); hy = 1/(npy+1); %gridsizes
bc    = 'som';
flag  = 1;

%% Test of the function klay.m
[x,y] = meshgrid(hx:hx:(1-hx),hy:hy:(1-hy));
k     = klay(x,y,kref);   
figure(1);
contourf(x,y,k);
%colormap(white)


%% Test of the function helmholtz2var.m 

% Wedge problem
kvar  = @(x,y) klay(x,y,kref);
eps  =  @(x,y) 0.5*(klay(x,y,kref).^2); 
zero  = @(x,y) 0*x;

[A]   = helmholtz2var(kvar,zero,npx,npy,bc);
[M]   = helmholtz2var(kvar,eps,npx,npy,bc);

% Random problem
%kvar  = @(x,y) krandn(x,y,kref);
%[A]  = helmholtz2var(kvar,npx,npy,bc);
%[M]   = shift_laplace2var(kvar,b1,b2,npx,npy,bc);


%% Solution of the Helmholtz equation and postprocessing the solution

if strcmp(bc,'dir')  
    [x,y]  = meshgrid(hx:hx:(1-hx),hy:hy:(1-hy));
    %right hand side: point source in the center of the domain
    b = zeros(length(A),1);
    b(npx*(npy-1)/2)= 1/(hx*hy);
    sol = A\b ;
    u   = reshape(sol',[npy,npx]);
   
else if strcmp(bc,'som')        
    [x,y] = meshgrid(0:hx:1,0:hy:1);
    
    %right hand side: point source
    b = zeros(length(A),1);
    index = ceil(npx+2)*(npy+1) + (npx/2);
    %index = ceil((npx+1)*(npy/2+1)+npy/2);
    
    b(index)= 1/(hx*hy);
   
    sol = A\b ;
    u = reshape(sol,[npx+2,npy+2]);
    
    end    
end
    
Reu  = real(u);  
Imu  = imag(u);

h = figure(2);
surf(y,x,Imu,'Edgecolor','none')
colormap()
set(h, 'linestyle', 'none');
