function [p] = gmres2D_csl_vs_dcsl_kvarying_coarse_inexact(kk,factoreps)

dim       = 2;
poweps    = 2;
m  = length(kk);
save_flag = 1;
p=1;

itercsl = zeros(m,1);
iterdef = zeros(m,3);

reflevs   = 1;   
restart   = [];
tol       = 1e-6;
maxit     = 200;
npcc      = 3;
par       = 0.7;


for i = 1:m
    
k    = kk(i);
bc   = 'som';
npf  = ceil(k^(3/2));
np   = npf-2;

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

[eqn1,~] = helmholtz2Dfem(node,elem,pdehelm,bdFlag,bdEdge);
[eqn2,~] = helmholtz2Dfem(node,elem,pdeSL,bdFlag,bdEdge);

[mg_mat,mg_split,restr,interp]= mg_setupfem_2D(npcc,numlev,pdeSL);

A    = eqn1.A;
Aeps = mg_mat{1};

n1 = length(A);
n2 = length(Aeps);

assert(n1==n2,'incorrect matrix size')

b       = ones(length(A),1);
restart = [];
u0      = zeros(length(A),1);

npre = 1; npos = 1; w  = 0.7; smo = 'wjac';
Aepsinv = @(x) Fcycle(mg_mat,mg_split,restr,interp,u0,x,npre,npos,w,smo,1);
mat1    = @(x) A*Aepsinv(x);

[~, ~, ~, iter1, ~] = gmres(mat1,  b, restart, tol, maxit);

itercsl(i)=iter1(2);
P = interp{1};
R = restr{1};
Ac = R*A*P;

setup.type = 'crout';
setup.droptol = 1e-1;
[Lc,Uc] = ilu(Ac,setup); 
length(Ac)

cgc =  @(x) P*(Uc\(Lc\(R*x)));
def =  @(x) Aepsinv(x-A*cgc(x))+cgc(x);
mat2 = @(x) A*def(x);
 
[~, ~, ~,iter2, ~] = gmres(mat2, b, restart, tol, maxit);

iterdef(i,1) = iter2(2);

Pc = P;
Rc = R;

for j = 2:3
    
    Pc = Pc*interp{j};
    Rc = restr{j}*Rc;
    Acc = R*A*P;
    [Lc,Uc] = lu(Acc);
    
    cgcc = @(x) Pccc*(Uc\(Lc\(Rcc*x)));
    def  = @(x) Aepsinv(x-A*cgc(x))+cgc(x);
    mat2 = @(x) A*def(x);
    [~, ~, ~,iter2, ~] = gmres(mat2, b, restart, tol, maxit);
    
    iterdef(i,j) = iter2(2);  
    length(Acc)
end

end


%% Postprocessing   %% 
epss   = num2str(10*factoreps);
table_data   = [kk',itercsl,iterdef(:,1),iterdef(:,2),iterdef(:,3)];
numrows      = size(table_data,1);
numcols      = size(table_data,2);
tableCaption = strcat('Number of GMRES iterations for the Helmholtz linear system preconditioned by the  CSL and the two-level method (TL) with different levels of coarsening');
tableLabel   = strcat('gmres_csl_vs_adef_coarse_eps_',epss);

%Data Format: .0f no decimals,
dataFormat = {'%.0f','%.0f','%.0f','%.0f','%.0f'};
header = ['\begin{tabular}','{',repmat(('c'),1,numcols),'}'];
table1  = {'\begin{table}[t]';'\centering';header;'\hline'};
row1   = {'$k$ & CSL & TL-1 & TL-2 & TL-3 \\ \hline'};
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
namefile = strcat('table_csl_vs_adef_coarse_facteps',epss,'.tex');
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