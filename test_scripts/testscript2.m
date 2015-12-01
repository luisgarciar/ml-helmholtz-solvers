%% Script for testing the function /matrixconstructors/helmholtz2d.m

%% Parameters
np   = 50;   %number of interior discretization points in 1D
npx  = np;
npy  = np;
hx   = 1/(npx+1); hy   = 1/(npy+1); %gridsizes
kmax = 10; kmin = 20;
flag = 1;  
dim  = 2;        %dimension    
bc   = 'dir';    %type of problem        
f    = @rhs2d;   %right hand side
k    = @(x,y)(krand(x,y,kmin,kmax)); %wavenumber (nonconstant)

if strcmp(bc,'dir')
    nv = npx*npy;    
else if strcmp(bc,'som')
    nv = (npx+2)*(npy+2);
    else
     error('invalid boundary conditions')
    end
end


%% Solution of the Helmholtz equation and postprocessing the solution
[A1, sol,b] = helmholtz2d(f,k,npx,npy,bc,flag);
[A2]        = helmholtz2(kref,npx,npy,bc);

b=ones(size(b));

x1=A1\b;
x2=A2\b;


if strcmp(bc,'dir')
    [x,y] = meshgrid(hx:hx:(1-hx),hy:hy:(1-hy));
       u  = reshape(sol',[npx,npy]);
   
else if strcmp(bc,'som')
    [x,y] = meshgrid(0:hx:1,0:hy:1);
        u = reshape(sol',[npx+2,npy+2]);
    end
    
end
    
Reu  = real(u);
Imu  = imag(u);

figure
%surf(x,y,Reu)
surf(x,y,feval(k,x,y))

