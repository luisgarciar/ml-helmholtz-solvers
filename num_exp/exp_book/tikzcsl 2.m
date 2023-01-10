function tikzcsl(k,ppw,b1,b2)
%% TIKZCSL Generates tikz files with plots of 
%  eigenvalues of the 1D Helmholtz
%  operator with Dirichlet boundary conditions
%  preconditioned by the shifted Laplacian and
%  the deflated CSL
%
%  Use:     tikzcsl(k,ppw,b1,b2)
%
%  Input:
%             k: Wavenumber
%             ppw: points per wavelength
%             (b1,b2): shift
%        
%  Output: 
%     csl_k_ppw_b1_b2.tex, dcsl_k_ppw_b1_b2.tex:   
%     files with the tikz figures (saved on current folder)     
%             
%  Author: Luis Garcia Ramos, 
%          Institut fur Mathematik, TU Berlin
%          Version 0.1, Feb 2016          
%%

clf;
[eCSL,eDCSL] = eigSL(k,ppw,b1,b2);

%setting file names for the .tex files
%filename format: csl_k_ppw_b1_b2.tex, dcsl_k_ppw_b1_b2.tex
wn     = num2str(k);  pts = num2str(ppw);
rshift = num2str(b1); ishift = num2str(b2);

plot(real(eCSL),imag(eCSL),'b.','MarkerSize',16);
axis equal
axis([-1 1 -1 1])
xlabel('Re(z)','FontSize',14) % x-axis label
ylabel('Im(z)','FontSize',14) % y-axis label
set(gca, 'FontSize', 12)
figure(1)
name1 = strcat('csl','_',wn,'_',pts,'.tex');

%save .tex file of tikz figure%
matlab2tikz('filename',name1,'standalone',true,'extraaxisoptions',...
           ['xlabel style={font=\Large},', ...
            'ylabel style={font=\Large},', ...
            'xtick distance={0.5},','ytick distance={0.5},']);


        
plot(real(eDCSL),imag(eDCSL),'b.','MarkerSize',16);
axis equal
axis([-1 1 -1 1])
xlabel('Re(z)','FontSize',14) % x-axis label
ylabel('Im(z)','FontSize',14) % y-axis label
set(gca, 'FontSize', 12)
name2 = strcat('dcsl','_',wn,'_',pts,'.tex');

%save .tex file of tikz figure%
matlab2tikz('filename',name2,'standalone',true,'extraaxisoptions',...
           ['xlabel style={font=\Large},', ...
            'ylabel style={font=\Large},', ...
            'xtick distance={0.5},','ytick distance={0.5},']);

        
        
       
end