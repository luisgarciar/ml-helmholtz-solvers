%% Script for testing the functions /rhs_k_constructors/klayer.m and krand.m

kmax  = 20; kmin = 30;
kref  = 40;
np    = 100;   %number of interior discretization points in 1D
npx   = np;
npy   = np;
hx    = 1/(npx+1); hy = 1/(npy+1); %gridsizes
b1    = 1;
b2    = 0.5;
bc    = 'som';
flag  = 1;

%% Test of the function klayer.m
[x,y] = meshgrid(hx:hx:(1-hx),hy:hy:(1-hy));

f     = @(x,y) rhs2d(x,y);
%kvar  = @(x,y) kref*ones(size(x));
kvar  = @(x,y) klayer(x,y,kref)
k     = klayer(x,y,kref);   
figure(1);
surf(x,y,k);

%break
%[x,y] = meshgrid(hx:hx:1-hx,hy:hy:1-hy);
% kk   = feval(kvar,x,y);
% kk   = reshape(kk',[nv,1]);


%% Test that the function helmholtz2.m works with variable wavenumber
[A1,sol,b] = helmholtz2d(f,kvar,npx,npy,bc,flag);
[A2]       = helmholtz2(kref,npx,npy,bc);

%x1=A1\b;
%x2=A2\b;


%% Solution of the Helmholtz equation and postprocessing the solution
[A,sol,b] = helmholtz2d(f,kvar,npx,npy,bc,flag);
%[M,~,~]   = shift_laplace2d(f,kvar,b1,b2,npx,npy,bc,flag);

if strcmp(bc,'dir')
    
    [x,y]  = meshgrid(hx:hx:(1-hx),hy:hy:(1-hy));
       b   = zeros(size(b));
       b(npx*(npy-1)/2)= 1; 
       sol = A1\b ;
       u   = reshape(sol',[npy,npx]);
   
else if strcmp(bc,'som')
        
    [x,y] = meshgrid(0:hx:1,0:hy:1);
     b = zeros(size(b));
     index = (npx+1)*(npy/2+1)+npy/2;
     b(index)= 1;
     sol = A1\b ;
     u = reshape(sol,[npy+2,npx+2]);
     %u = reshape(b,[npy+2,npx+2]);

    end    
end
    
Reu  = real(u);
Imu  = imag(u);

figure(2)
surf(x,y,Reu)