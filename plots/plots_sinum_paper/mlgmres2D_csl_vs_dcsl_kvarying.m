function [p] = mlgmres2D_csl_vs_dcsl_kvarying(kk,factoreps)

dim       = 2;
poweps    = 2;
m         = length(kk);
save_flag = 1;
p         = 1;

itercsl   =  zeros(m,1);
iterdef   =  zeros(m,2);
timecsl   =  zeros(m,1);
timedef   =  zeros(m,2);

tol       =  1e-6;
maxit     =  200;
npcc      =  1;
par       =  0.6;
bc        = 'som';


for i = 1:m
    k    = kk(i);
    [npf,numlev] = fem_npc_to_npf(npcc,k,par);  %number of points in finest grid (1D)    
    h = 1/(npcc+1);
    [node,elem] = squaremesh([0 1 0 1],h);
    
    %refining the mesh numlev times
    for j = 1:numlev-1
        [node,elem] = uniformrefine(node,elem);
    end
    
    %set boundary conditions
    [bdNode,bdEdge,isBdNode] = findboundary(elem);
    bdFlag = setboundary(node,elem,'ABC');
    
    pdehelm    = helmholtz2Dconstantwndata(k,0,1);
    pdeSL      = helmholtz2Dconstantwndata(k,factoreps,poweps);
    option.tol = 1e-8;
    
    %building the matrices
    [eqn1,~] = helmholtz2Dfem(node,elem,pdehelm,bdFlag,bdEdge);
    [eqn2,~] = helmholtz2Dfem(node,elem,pdeSL,bdFlag,bdEdge);
    
    %setting up the multilevel structure
    [mg_matHelm,mg_splitHelm,restr,interp]= mg_setupfem_2D(npcc,numlev,pdehelm);
    [mg_matCSL,mg_splitCSL,restrCSL,interpCSL]= mg_setupfem_2D(npcc,numlev,pdeSL);
    
    A      = mg_matHelm{1};
    Aeps   = mg_matCSL{1};
    
    n1 = length(A);
    n2 = length(Aeps);
    
    assert(n1==n2,'incorrect matrix size')
    
    b       = ones(length(A),1);
    restart = [];
    u0      = zeros(length(A),1);
    
    npre = 1; npos = 1; w  = 0.7; smo = 'wjac';
    Aepsinv = @(x) Fcycle(mg_matCSL,mg_splitCSL,restrCSL,interpCSL,u0,x,npre,npos,w,smo,1);
    AP  = @(x) A*Aepsinv(x);
    tol = 1e-6;
    
    
    factorepss = num2str(factoreps);
    ks = num2str(kk(i));
    
    %% CSL Run
    msg = strcat('Beginning GMRES-CSL run for Helmholtz problem with k','=',ks,', eps','=',factorepss,'*k^2');
    disp(msg);
    disp(strcat('Size of problem is',': ',num2str(length(A))));
    
    tic
    [~, ~, ~, iter1, ~] = gmres(AP,b,restart,tol,maxit);
    timecsl(i,1) = toc;
    
    msg = strcat('Finished GMRES-CSL run');
    disp(msg);
    
    itercsl(i)=iter1(2);
    
    %% ML Run with MK(8,4,2,1)
    maxiter = ones(length(mg_matHelm),1);
    maxiter(1:5)=[20,8,2,2,1]';
    x0 = zeros(n1,1);
    
    msg = strcat('Beginning GMRES-TL run for Helmholtz problem with k', '=',ks,', eps', ' = ', factorepss,'*k^2');
    tic
    [x2,~,~,iter2] = mlfgmres(b,x0,mg_matHelm,mg_matCSL,mg_splitCSL,restr,interp,maxiter,tol);
    timeml1 = toc;
    msg = strcat('Finished GMRES-TL run');
    iterdef(i,1) = iter2;
    timeml(i,1)=timeml1;
    
    %% ML Run with MK(6,4,2,1)

    maxiter(1:5)=[20,6,2,2,1]';
    
    msg = strcat('Beginning GMRES-TL run for Helmholtz problem with k', '=',ks,', eps', ' = ', factorepss,'*k^2');
    disp(msg)
    tic
    [x2,~,~,iter2] = mlfgmres(b,x0,mg_matHelm,mg_matCSL,mg_splitCSL,restr,interp,maxiter,tol);
    timeml2 = toc;
    msg = strcat('Finished GMRES-TL run');
    disp(msg);
    
    iterdef(i,2) = iter2;
    timeml(i,2)=timeml2;

end


%% Postprocessing   %%
epss = num2str(10*factoreps);
table_data   = [kk',itercsl,timecsl,iterdef(:,1),timeml(:,1),iterdef(:,2),timeml(:,2)];
numrows      = size(table_data,1);
numcols      = size(table_data,2);
tableCaption = strcat('Number of GMRES iterations and total computation time (in seconds) for the Helmholtz linear system preconditioned by the  CSL and the multilevel method (ML)');
tableLabel   = strcat('mlgmres_csl_vs_adef_coarse_eps_',epss);

%Data Format: .0f no decimals,
dataFormat = {'%.0f','%.0f','%g','%.0f','%g','%.0f','%g'};
header = ['\begin{tabular}','{',repmat(('c'),1,numcols),'}'];
table1  = {'\begin{table}[t]';'\centering';header;'\hline'};
row1   = {'$k$ & CSL & CSL-t & ML(8,4,2) & ML(8,4,2)-t& ML(6,4,2) & ML(6,4,2)-t \\ \hline'};
table1 = [table1;row1'];

for i=1:numrows
    for j =1:numcols
        dataValue = num2str(table_data(i,j),dataFormat{j});
        if j==1
            rowStr = dataValue;
        else
            rowStr = [rowStr,' & ',dataValue];
        end
    end
    table1(end+1) = {[rowStr,' \\']};
end

footer = {'\end{tabular}';['\caption{',tableCaption,'}']; ...
    ['\label{table:',tableLabel,'}'];'\end{table}'};
table1 = [table1;'\hline';footer];

%Save the table to a file in folder ../../../tex_files/tables
namefile = strcat('table_mlgmres_csl_vs_ml_facteps',epss,'.tex');
%currentpath = mfilename('fullpath');
currentpath = fileparts(mfilename('fullpath'));

f = fullfile(currentpath,'/plots/tex/tables',namefile);

if (save_flag == 1)
    fid = fopen(f,'w');
    [nrows,ncols] = size(table1);
    for row = 1:nrows
        fprintf(fid,'%s\n',table1{row,:});
    end
    fclose(fid);
end