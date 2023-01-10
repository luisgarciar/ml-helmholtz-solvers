%% Experiments with the CSL-preconditioned Helmholtz equation.
% We compare timing and the iteration numbers for solving
%    A M^{-1}x = b
% with gmres, for different implementations of A*M^{-1}*x:
%  - eps = k
%  - eps = k^2

clear all;
save_flag = 1;  % save_flag=1: save plots and table, =0 do not save.
% Setup parameters
% Setup list of wavenumbers and shifts
%wavenum   = [20 40 60] ;
%wavenum   = [10 20 30];
%wavenum = 0;
wavenum = [10 20 40 60 80 100 120];
poweps = 1;

%wavenum = 50; %% run this when testing changes in the code
%number of interior points in coarsest grid in one dim
npc = 1;
bc = 'som';
ppw = 12;   

linetyp = {'-','-.','--'};
color  = [0 0 0; 0.5 0 0.5; 0 0 1; 0 0.5 0; 1 0 0];

%Parameters for GMRES
restart   = [];
tol       = 1e-8;
maxit     = 150;

pollution = 'no';

%memory allocation for iteration counts
iter_num1D = zeros(length(wavenum),length(poweps));
iter_num2D = zeros(length(wavenum),length(poweps));


for kk = 1:length(wavenum)
    %wavenumber
    k = wavenum(kk);
    eps = k^poweps;
    
    %number of points in finest grid and number of levels
    %choosing the number of points
    
    if strcmp(pollution,'no')
        npf = ceil(k^(3/2));
    end
    
    if (mod(npf+1,2)==1)  %set an even number of interior points in 1D
        npf = npf+1;
    end
    npc = (npf-1)/2;
    
    %1D
    A    = helmholtzfem(k,npf,0,bc);           %Helmholtz matrix
    Aeps = helmholtzfem(k,npf,eps,bc);         %Shifted Laplace matrix
    b    = ones(length(A),1);
    
    [L, U] = lu(Aeps); % LU = M
    mat  = @(x) A*(U\(L\x));
    
    [~, ~, ~, iter, ~] = gmres(mat, b, restart, tol, maxit);
    iter_num1D(kk,1) = iter(2);
    
    %2D
    dim  = 2;
    A    = helmholtz2(k,0,npf,npf,bc);
    Aeps = helmholtz2(k,eps,npf,npf,bc);
    op_type = 'gal';
    numlev = 2;
    %[mg_mat,mg_split,restrict,interp] = mg_setup(k,eps,op_type,npc,numlev,bc,dim);
    
    b  = ones(length(A),1);
    x0 = zeros(size(b));
    
    if k<=60
        [L, U] = lu(Aeps); % LU = M
        mat    = @(x) A*(U\(L\x));        
        [~, ~, ~, iter, ~] = gmres(mat, b, restart, tol, maxit);
        iter_num2D(kk,1)   = iter(2);
        
    else
        iter_num2D(kk,1) = 0;
    end
end

%% Plot iteration numbers and save as tikz (.tex) and .eps
figure(2)

%plot iteration counts 
%preconditioned problem 1D
plot(wavenum, iter_num1D(:,1),'k-','Linewidth',5);
hold on

iter_num2D60 = iter_num2D(wavenum<=60);
%preconditioned problem_2D
plot(wavenum(wavenum<=60), iter_num2D60(:,1),'b-','Linewidth',5);

hold off

%set(gca,'Xtick',[0 0.5 1],'FontSize',30);
%set(gca,'Ytick',[-0.5 0 0.5],'FontSize',30);
FS = 20; % font size

ylabel('Number of GMRES iterations','fontsize',FS)
xlabel('Wavenumber','fontsize',FS)
%legend('No Prec.', ... 
legend('1D', ...
       '2D',...
       'Location','NorthWest','Fontsize',FS)
%title(['1D Helmholtz with CSL-preconditioner and poly. accel., deg=',...
%num2str(deg)])
%set(gca,'LooseInset',get(gca,'TightInset'))
%set(gca,'FontSize',FS)

axis([10 100 0 10])
set(gca,'Xtick',10:20:100,'FontSize',FS);
set(gca,'Ytick',0:2:10,'FontSize',FS);