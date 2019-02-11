%Plots of FOVs for SINUM paper


%% Plots with eps = k^2

%% Plots with N = ceil(k^2)
%1D Case
opts.dim = 1;
opts.poweps    = 2;
opts.factoreps = 1;
opts.bc = 'som';
opts.fvpts = 64;
opts.disc = 'q';

kmult  =  [5 20 50 80 100];

%CSL
opts.prec = 'csl';
[fovAAeps] = computefovAhat(kmult,opts);
plotfovAhat(fovAAeps,kmult,opts);
close all

%Deflated CSL
opts.prec = 'adef';
[fovAB] = computefovAhat(kmult,opts);
plotfovAhat(fovAB,kmult,opts);
close all



%% Plots with  N = ceil(k^1.5)
%1D Case
opts.dim = 1;
opts.poweps    = 2;
opts.factoreps = 1;
opts.bc = 'som';
opts.fvpts = 64;
opts.disc = 'pf';

kmult  =  [5 20 50 80 100];

%CSL
opts.prec = 'csl';
[fovAAeps] = computefovAhat(kmult,opts);
plotfovAhat(fovAAeps,kmult,opts);
close all

%Deflated CSL
opts.prec = 'adef';
[fovAB] = computefovAhat(kmult,opts);
plotfovAhat(fovAB,kmult,opts);
close all




%% Plots with eps = 5k^2

%% Plots with N = ceil(k^2)
%1D Case
opts.dim = 1;
opts.poweps    = 2;
opts.factoreps = 5;
opts.bc = 'som';
opts.fvpts = 64;
opts.disc = 'q';

kmult  =  [5 20 50 80 100];

%CSL
opts.prec = 'csl';
[fovAAeps] = computefovAhat(kmult,opts);
plotfovAhat(fovAAeps,kmult,opts);
close all

%Deflated CSL
opts.prec = 'adef';
[fovAB] = computefovAhat(kmult,opts);
plotfovAhat(fovAB,kmult,opts);
close all



%% Plots with  N = ceil(k^1.5)

%1D Case
opts.dim = 1;
opts.poweps    = 2;
opts.factoreps = 5;
opts.bc = 'som';
opts.fvpts = 64;
opts.disc = 'pf';

kmult  =  [5 20 50 80 100];

%CSL
opts.prec = 'csl';
[fovAAeps] = computefovAhat(kmult,opts);
plotfovAhat(fovAAeps,kmult,opts);
close all

%Deflated CSL
opts.prec = 'adef';
[fovAB] = computefovAhat(kmult,opts);
plotfovAhat(fovAB,kmult,opts);
close all




%% Plots with eps = 10k^2

%% Plots with N = ceil(k^2)
%1D Case
opts.dim = 1;
opts.poweps    = 2;
opts.factoreps = 10;
opts.bc = 'som';
opts.fvpts = 64;
opts.disc = 'q';

kmult  =  [5 20 50 80 100];

%CSL
opts.prec = 'csl';
[fovAAeps] = computefovAhat(kmult,opts);
plotfovAhat(fovAAeps,kmult,opts);
close all

%Deflated CSL
opts.prec = 'adef';
[fovAB] = computefovAhat(kmult,opts);
plotfovAhat(fovAB,kmult,opts);
close all



%% Plots with  N = ceil(k^1.5)

%1D Case
opts.dim = 1;
opts.poweps    = 2;
opts.factoreps = 10;
opts.bc = 'som';
opts.fvpts = 64;
opts.disc = 'pf';

kmult  =  [5 20 50 80 100];

%CSL
opts.prec = 'csl';
[fovAAeps] = computefovAhat(kmult,opts);
plotfovAhat(fovAAeps,kmult,opts);
close all

%Deflated CSL
opts.prec = 'adef';
[fovAB] = computefovAhat(kmult,opts);
plotfovAhat(fovAB,kmult,opts);
close all


%% 
% %% Plot and save  FOV CSL
% close all
% for i=1:length(kk)
%     k   = kk(i);   
%     eps = factoreps*k^poweps;
% 
%     label = strcat('$k = ', num2str(kmult(i)),' \pi$');
%     fovplotCSL(i) =  plot(real(fovAAeps(:,i)),imag(fovAAeps(:,i)),'Color',color1(i,:),...
%                       'LineWidth',2,'linestyle',linetyp{i},...
%                       'DisplayName',label);
%        
%     hold on
%     plot(0,0,'k.','Markersize',10,'LineWidth',1);
%     plot(1,0,'b.','Markersize',10,'LineWidth',1);
%     axis equal
%     axis([-0.2 1.2 -0.7 0.7]);
%     xlabel('Re(z)','FontSize',16,'Interpreter','latex');
%     ylabel('Im(z)','FontSize',16,'Interpreter','latex');
%     h=gca;
%    
%     set(gca,'Xtick',[0 0.5 1],'FontSize',16);
%     set(gca,'Ytick',[-0.5 0 0.5],'FontSize',16);
%     set(gca,'TickLabelInterpreter', 'tex');
% 
% end
% 
% L=legend(fovplotCSL);
% set(L,'Interpreter','latex','FontSize',16);
% 
% 
% kmin  = num2str(min(kmult));
% kmax  = num2str(max(kmult));
% pts   = num2str(ppw);
% polexp = num2str(10*pollexp);
% powershift  = num2str(poweps);
% factorshift = num2str(10*factoreps);
% 
% %Tikz Axis formatting
% x = xlabel('$\mathrm{Re}(z)$');
% %x-axis label
% set(x,'Interpreter','latex','fontsize',16)
% y=ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',16); % x-axis label
% set(y,'Interpreter','latex','fontsize',16)
% fig = figure(1);
% 
% nametex = strcat('1d_fov_csl_kmin',kmin,'_kmax',kmax,'_ppw',pts, ...
%     '_pshift_',powershift,'_10timesfacshiftfshift_',factorshift,'.tex');
% 
% if strcmp(pollution,'no')
%     nametex = strcat('1d_fov_csl_kmin',kmin,'_kmax',kmax,'_nopoll10exp_',polexp, ...
%         '_powshift_',powershift,'_10timesfacshift_',factorshift,'.tex');
% end
% 
% %path for saving the files
% currentpath = pwd;
% filetex = fullfile(currentpath,'plots','tex',nametex);
% 
% %save .tex file of tikz figure%
% matlab2tikz('filename',filetex,'standalone',true,...
%              'interpretTickLabelsAsTex',true,...
%              'extraaxisoptions',['xlabel style={font=\Large},'...
%              'ylabel style={font=\Large},',...
%              'legend style={font=\Large},',...
%              'ticklabel style={font=\Large}']);
% 
% nameeps = strcat('1d_fov_csl_kmin',kmin,'_kmax',kmax,'_ppw',pts, ...
%     '_pshift_',powershift,'_10timesfacshift_',factorshift,'.eps');
% %nameeps = 'something';
% 
% plot_file_eps = fullfile(currentpath,'plots','eps',nameeps);
% print(fig,'-depsc',plot_file_eps)         
%          
% %% Plot and save FOV Deflated CSL
% close all
%         
% for i=1:length(kk)
%     k   = kk(i);   
%     eps = factoreps*k^poweps;
% 
%     label = strcat('$k = ', num2str(kmult(i)),' \pi$');
%     fovplotdefCSL(i) =  plot(real(fovAB(:,i)),imag(fovAB(:,i)),'Color',color1(i,:),...
%                       'LineWidth',2,'linestyle',linetyp{i},...
%                       'DisplayName',label);
%        
%     hold on
%     plot(0,0,'k.','Markersize',10,'LineWidth',1);
%     plot(1,0,'b.','Markersize',10,'LineWidth',1);
%     axis equal
%     axis([-0.2 1.2 -0.7 0.7]);
%     xlabel('Re(z)','FontSize',16,'Interpreter','latex');
%     ylabel('Im(z)','FontSize',16,'Interpreter','latex');
%     h=gca;
%    
%     set(gca,'Xtick',[0 0.5 1],'FontSize',16);
%     set(gca,'Ytick',[-0.5 0 0.5],'FontSize',16);
%     set(gca,'TickLabelInterpreter', 'tex');
% 
% end
% 
% L=legend(fovplotdefCSL);
% set(L,'Interpreter','latex','FontSize',16);
% 
% 
% kmin  = num2str(min(kmult));
% kmax  = num2str(max(kmult));
% pts   = num2str(ppw);
% polexp = num2str(10*pollexp);
% powershift  = num2str(poweps);
% factorshift = num2str(10*factoreps);
% 
% %Tikz Axis formatting
% x = xlabel('$\mathrm{Re}(z)$');
% %x-axis label
% set(x,'Interpreter','latex','fontsize',16)
% y=ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',16); % x-axis label
% set(y,'Interpreter','latex','fontsize',16)
% fig = figure(1);
% 
% nametex = strcat('1d_fov_defcsl_kmin',kmin,'_kmax',kmax,'_ppw',pts, ...
%     '_pshift_',powershift,'_10timesfacshiftfshift_',factorshift,'.tex');
% 
% if strcmp(pollution,'no')
%     nametex = strcat('1d_fov_defcsl_kmin',kmin,'_kmax',kmax,'_nopoll10exp_',polexp, ...
%         '_powshift_',powershift,'_10timesfacshift_',factorshift,'.tex');
% end
% 
% %path for saving the files
% currentpath = pwd;
% filetex = fullfile(currentpath,'plots','tex',nametex);
% 
% %save .tex file of tikz figure%
% matlab2tikz('filename',filetex,'standalone',true,...
%              'interpretTickLabelsAsTex',true,...
%              'extraaxisoptions',['xlabel style={font=\Large},'...
%              'ylabel style={font=\Large},',...
%              'legend style={font=\Large},',...
%              'ticklabel style={font=\Large}']);
% 
% nameeps = strcat('1d_fov_defcsl_kmin',kmin,'_kmax',kmax,'_ppw',pts, ...
%     '_pshift_',powershift,'_10timesfacshift_',factorshift,'.eps');
% %nameeps = 'something';
% 
% plot_file_eps = fullfile(currentpath,'plots','eps',nameeps);
% print(fig,'-depsc',plot_file_eps) 
%   

