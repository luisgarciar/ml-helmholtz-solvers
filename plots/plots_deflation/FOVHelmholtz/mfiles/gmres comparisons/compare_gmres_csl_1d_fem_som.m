%% Experiments with the CSL-preconditioned Helmholtz equation.
% We compare timing and the iteration numbers for solving
%    A M^{-1}x = b
% with gmres, for different implementations of A*M^{-1}*x:
%  - no prec, M=I
%  - eps = k
%  - eps = k^2
%
clear all;
save_flag = 1;  % save_flag=1: save plots and table, =0 do not save.

% Setup parameters
% Setup list of wavenumbers and shifts
wavenum   = 20:20:100 ;
poweps    = [1 1.5 2];

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
maxit     = 100;

pollution = 'no';

%memory allocation for iteration counts
iter_num = zeros(length(wavenum),length(poweps));


for kk = 1:length(wavenum)
    %wavenumber
    k = wavenum(kk);
    
    for j=1:length(poweps)     
        eps = k^poweps(j);
        
        %number of points in finest grid and number of levels
        %choosing the number of points
        npf = ceil(ppw*k/(2*pi))-1;
        
        if strcmp(pollution,'no')
            npf = ceil(k^(3/2));
        end
        
        if (mod(npf+1,2)==1)  %set an even number of interior points in 1D
            npf = npf+1;
        end
        npc = (npf-1)/2;
        
        %1D matrices
        A    = helmholtzfem(k,npf,0,bc);           %Helmholtz matrix
        Aeps = helmholtzfem(k,npf,eps,bc);         %Shifted Laplace matrix
        
        %2D
        %dim  = 2;
        %A    = helmholtz2(k,0,npf,npf,bc);
        %Aeps = helmholtz2(k,eps,npf,npf,bc);
        %op_type = 'gal';
        %[mg_mat,mg_split,restrict,interp] = mg_setup(k,eps,op_type,npc,numlev,bc,dim);
        
        b  = ones(length(A),1);
        x0 = zeros(size(b));
                
        if (poweps(j)~=0)
        [L, U] = lu(Aeps); % LU = M
        mat  = @(x) A*(U\(L\x));
        else
        mat  = @(x) A*x;
        end 
        
        [~, ~, ~, iter, ~] = gmres(mat, b, restart, tol, maxit);
        iter_num(kk,j) = iter(2);
        
    end% of going through different poweps
    
end % of going through different wavenumbers


%% Plot iteration numbers and save as tikz (.tex) and .eps
figure(2)

%plot iteration counts 
%preconditioned problem k=1
plot(wavenum, iter_num(:,1),'Color',color(1,:),...
    'linestyle',linetyp{1},'Linewidth',5);
hold on

%preconditioned problem
for j=2:length(poweps)
plot(wavenum, iter_num(:,j), 'Color',color(j,:),...
 'Linestyle',linetyp{j},'Linewidth',4);

end

%set(gca,'Xtick',[0 0.5 1],'FontSize',30);
%set(gca,'Ytick',[-0.5 0 0.5],'FontSize',30);


hold off
ylabel('Number of GMRES iterations','fontsize',25)
xlabel('Wavenumber','fontsize',25)
%legend('No Prec.', ... 
legend('\epsilon=k', ...
     '\epsilon=k^{3/2}',...
     '\epsilon=k^2',... 
     'Location','NorthEast')
%title(['1D Helmholtz with CSL-preconditioner and poly. accel., deg=',...
%num2str(deg)])
FS = 25; % font size
set(gca,'LooseInset',get(gca,'TightInset'))
set(gca,'FontSize',FS)

plot_iter_tex = strcat('comp_iter_vs_k_eps_1D_som.tex');

%get path of current .m file
currentpath  = mfilename('fullpath');
%generate file name, for saving in tex_files/figures/new_exp
plot_file_tex = fullfile(currentpath,plot_iter_tex);

if ( save_flag == 1 )
    %Save as tikz figure in .tex file
    matlab2tikz('filename',plot_iter_tex,'standalone',true,'extraaxisoptions',...
                ['xlabel style={font=\LARGE},', ...
                 'ylabel style={font=\LARGE},']);
    
    %Save the plot as .eps
   % plot_iter_eps = strcat('comp_iter_vs_k_eps_1D_som.eps');
   % plot_file_eps = fullfile(path,plot_iter_eps);
   % print('-depsc2',plot_file_eps)
end


