%% test of linear interpolation
 clear all;
 npc = 2; dim=2;  bc = 'som';
 npcx  = npc; npcy = npc;
 npccx = npcx+2; npccy= npcy+2;
 npffx = 2*npccx-1; npffy = 2*npccy-1; 
 M = npffx*npffy; N = npccx*npccy;

 Z = sparse(npffx,npccx); %z(1:3,1) = [1;2;1];
 u = 0.5*ones(npccx,1);
 odd   = (1:2:npffx);  %indices of coarse grid points
 even  = (2:2:npffx-1);
 
 Z1 = spdiags([u u],[1 0],npccx-1,npccx);
 Z  = sparse(npffx,npccx);
 Z(odd,:) = speye(npccx);
 Z(even,:)= Z1;
 Z2d = kron(Z,Z); 
 
 %Z   = gallery('circul',z)';
 %Z   = 0.5*sparse(Z(:,1:2:(npccx)));
 
 
 %% Other stuff
 
%If fine grid point coincides with coarse grid point
%transfer the function value unchanged from coarse to fine
i  = [1:2:npffx]; ic = ceil(i/2);
j  = [1:2:npffy]; jc = ceil(j/2);
Z  = sparse(npffx,npffy); Z(i,j)=1;    Z  = reshape(Z',npffx*npffy,1);
Z2 = sparse(npccx,npccy); Z2(ic,jc)=1; Z2 = reshape(Z2',npccx*npccy,1);

%R: (fine grid) indices of the rows to be changed
%C: C(i) is the column (coarse grid index) to be filled in row R(i)
%R*numcolumns + C has the index of the entry to be changed
R   = find(Z);
C   = (1:1:N)';
%Matrix entries to be modified
ind     = (R-1)*N+C; 
Interp1 = sparse(M*N,1); Interp1(ind)=1; 
Interp1 = full(reshape(Interp1',N,M))'; 
 
%If a fine grid point has only south and north coarse grid 
%neighbors, interpolate from function values at neighbor points
%using weights 0.5
i    = 1:2:npffx;     
j    = 2:2:(npffy-1);

%R:Row Indices of fine grid points
Z = sparse(npffx,npffy); Z(i,j)=1; Z=reshape(Z,npffx*npffy,1)';
R  = find(Z); 
%CS,CN: Indices of south and north coarse grid neighbors
CS =  1:1:(npccx*(npccy-1));
CN = (npccx+1:1:npccx*npccy);

indS = (R-1)*N+CS;
indN = (R-1)*N+CN;

Interp2 = sparse(M*N,1); Interp2(indS)=0.5; Interp2(indN)=0.5;
Interp2 = reshape(Interp2',N,M)'; 

%If a fine grid point has only east and west coarse grid 
%neighbors, interpolate from function values at neighbor points
%using weight 0.5
i = 2:2:npffx-1; ic = ceil(i/2);
j = 1:2:npffy;   jc = ceil(j/2);

%R: Row Indices of fine grid points with only east and west coarse neighbors
Z = sparse(npffx,npffy); Z(i,j)=1; Z=reshape(Z,npffx*npffy,1)';
R  = find(Z); 

%CS,CN: Column Indices of west and east coarse grid neighbors
iw = 1:1:(npccx-1);
jw = 1:1:npccy;
ZSW = sparse(npccx,npccy); ZSW(iw,jw)=1; ZSW=reshape(ZSW,npccx*npccy,1)';
CW = find(ZSW);
CE = CW+1;

indW = (R-1)*N+CW;
indE = (R-1)*N+CE;

Interp3 = sparse(M*N,1); Interp3(indW)=0.5; Interp3(indE)=0.5;
Interp3 = reshape(Interp3',N,M)'; 

%If a fine grid point has four coarse grid neighbors
% interpolate from function values at neighbor points
%using weights 0.25

i = 2:2:npffx-1; 
j = 2:2:npffy-1;

Z  = sparse(npffx,npffy); Z(i,j)=1; Z=reshape(Z,npffx*npffy,1)';
R = find(Z); 

%South-West neighbors
iSW = 1:1:(npccx-1);
jSW = 1:1:(npccy-1);
[I,J] = meshgrid(iSW,jSW);
ZSW = sparse(npccx,npccy); ZSW(I,J)=1; 
ZSW=reshape(ZSW',npccx*npccy,1)';
CSW = find(ZSW);

%Other Neighbors
CSE = CSW+1;
CNW = CSW+npccx; CNE = CNW+1;

%Indices of matrix entries to be modified
ind = ([R R R R]-1)*N+[CSW CSE CNW CNE];

Interp4 = sparse(M*N,1); Interp4(ind)=0.25;
Interp4 = reshape(Interp4',N,M)'; 

%Interpolation matrix
Interp  = Interp1+Interp2+Interp3+Interp4;
Interpf = lininterpol(npc,dim,bc);



%[I J] = meshgrid(i,j);

%grid = [npffx npffy];
%R1   = reshape(sub2ind(grid,J,I),1,[]);


%  for j=1:npffx
%      for i=1:npffy
%          indf = i+(j-1)*npffx;          
%         if (mod(i,2)==1 && mod(j,2)==1)
%            %if fine grid point coincides with coarse grid point
%            %transfer the function value unchanged
%             %sequential index of point
%            ii = intdiv(i,2); jj=intdiv(j,2); 
%            indc = ii+(jj-1)*npccx;
%            Z(indf,indc)=1;
% 
%         elseif (mod(i,2)==0 && mod(j,2)==1)
%             %interpolate using east and west neighbors from
%             %coarse grid
%             iiwest = intdiv(i,2); jj=intdiv(j,2);
%             %index of west neighbor in coarse grid
%             indwestc = iiwest+(jj-1)*npccx; 
%             ind = [indwestc indwestc+1];
%             assert(max(ind)<=N,'index too large'); 
%             Z(indf,ind) = 0.5; 
% 
%         elseif(mod(i,2)==1 && mod(j,2)==0)
%                %interpolate using south and north coarse grid
%                %neighbors
%                iisouth   = intdiv(i,2); jjsouth=intdiv(j,2);
%                indsouthc = iisouth+(jjsouth-1)*npccx;
%                ind = [indsouthc indsouthc+npccx];
%                assert(max(ind)<=N,'index too large'); 
%                Z(indf,ind) = 0.5;
% 
%         elseif(mod(i,2)==0 && mod(j,2)==0)
%             %interpolate using SW,SE,NW,NE coarse grid
%             %neighbors
%             iiSW = intdiv(i,2); jjSW = intdiv(j,2);
%             iiSE = iiSW+1;      jjSE = jjSW;
%             iiNW = iiSW;        jjNW = jjSW+1;
%             iiNE = iiSE;        jjNE = jjSW+1;
%             ii = [iiSW iiSE iiNW iiNE];
%             jj = [jjSW jjSE jjNW jjNE];
%             ind = ii+(jj-1)*npccx;
%             assert(max(ind)<=N,'index is too large');
%             Z(indf,ind) = 0.25;
%         end
%      end
%  end
%    



%OLD CODE
%                %If fine grid point coincides with coarse grid point
%                %transfer the function value unchanged from coarse to fine
%                i  = [1:2:npffx]; ic = ceil(i/2);
%                j  = [1:2:npffy]; jc = ceil(j/2);
%                Z1  = sparse(npffx,npffy); Z1(i,j)=1;  Z1  = reshape(Z1',npffx*npffy,1);
%                Z2  = sparse(npccx,npccy); Z2(ic,jc)=1; Z2 = reshape(Z2',npccx*npccy,1);
% 
%                %R: (fine grid) indices of the rows to be changed
%                %C: C(i) is the column (coarse grid index) to be filled in row R(i)
%                %R*numcolumns + C has the index of the entry to be changed
%                R   = find(Z1);
%                C   = (1:1:N)';
%                %Matrix entries to be modified
%                ind     = (R-1)*N+C; 
%                Interp1 = sparse(M*N,1); Interp1(ind)=1; 
%                Interp1 = reshape(Interp1',N,M)'; 
%                
%                %%If a fine grid point has only south and north coarse grid 
%                %neighbors, interpolate from function values at neighbor points
%                %using weights 0.5
%                i  = 1:2:npffx;     
%                j  = 2:2:(npffy-1);
% 
%                %R:Row Indices of fine grid points
%                Z1 = sparse(npffx,npffy); Z1(i,j)=1; Z1=reshape(Z1,npffx*npffy,1)';
%                R  = find(Z1); 
%                %CS,CN: Indices of south and north coarse grid neighbors
%                CS =  1:1:(npccx*(npccy-1));
%                CN = (npccx+1:1:npccx*npccy);              
%                indS = (R-1)*N+CS;
%                indN = (R-1)*N+CN;         
%                Interp2 = sparse(M*N,1); Interp2(indS)=0.5; Interp2(indN)=0.5;
%                Interp2 = reshape(Interp2',N,M)';
%                 
%                %If a fine grid point has only east and west coarse grid 
%                %neighbors, interpolate from function values at neighbor points
%                %using weight 0.5
%                i = 2:2:npffx-1; 
%                j = 1:2:npffy;  
% 
%                %R: Row Indices of fine grid points with only east and west coarse neighbors
%                Z1 = sparse(npffx,npffy); Z1(i,j)=1; Z1=reshape(Z1,npffx*npffy,1)';
%                R  = find(Z1); 
% 
%                %CS,CN: Column Indices of west and east coarse grid neighbors
%                iw = 1:1:(npccx-1);
%                jw = 1:1:npccy;
%                ZSW = sparse(npccx,npccy); ZSW(iw,jw)=1; ZSW=reshape(ZSW,npccx*npccy,1)';
%                CW = find(ZSW);
%                CE = CW+1;
% 
%                indW = (R-1)*N+CW;
%                indE = (R-1)*N+CE;
% 
%                Interp3 = sparse(M*N,1); Interp3(indW)=0.5; Interp3(indE)=0.5;
%                Interp3 = reshape(Interp3',N,M)'; 
%                
%                %If a fine grid point has four coarse grid neighbors
%                % interpolate from function values at neighbor points
%                %using weights 0.25
% 
%                i = 2:2:npffx-1; 
%                j = 2:2:npffy-1;
% 
%                Z1  = sparse(npffx,npffy); Z1(i,j)=1; Z1=reshape(Z1,npffx*npffy,1)';
%                R = find(Z1); 
% 
%                %South-West neighbors
%                iSW = 1:1:(npccx-1);
%                jSW = 1:1:(npccy-1);
%                [I,J] = meshgrid(iSW,jSW);
%                ZSW = sparse(npccx,npccy); ZSW(I,J)=1; 
%                ZSW=reshape(ZSW',npccx*npccy,1)';
%                CSW = find(ZSW);
% 
%                %Other Neighbors
%                CSE = CSW+1;
%                CNW = CSW+npccx; CNE = CNW+1;
% 
%                %Indices of matrix entries to be modified
%                ind = ([R R R R]-1)*N+[CSW CSE CNW CNE];
% 
%                Interp4 = sparse(M*N,1); Interp4(ind)=0.25;
%                 Interp4 = reshape(Interp4',N,M)'; 
% 
%                 %Interpolation matrix
%                 Z  = Interp1+Interp2+Interp3+Interp4;          
%                
% 
%                %We fill the matrix by rows
%                for j=1:npffy
%                    for i=1:npffx
%                        %(i,j): lexicographic index of point (hx*(i-1),hy*(j-1))
%                        %in fine grid
%                        indf = i+(j-1)*npffx; %sequential index of point
%                        if (mod(i,2)==1 && mod(j,2)==1)
%                        %if fine grid point coincides with coarse grid point
%                        %transfer the function value unchanged
%                            ii = ceil(i/2); jj = ceil(j/2);
%                            indc = ii+(jj-1)*npccx;
%                            assert(indc<=N,'index too large')
%                            Z(indf,indc) = 1;
%                            
%                        elseif (mod(i,2)==0 && mod(j,2)==1)
%                            %interpolate using east and west coarse grid
%                            %neighbors 
%                            iiwest = ceil(i/2); jjwest=ceil(j/2);
%                            %index of west neighbor in coarse grid
%                            indwestc = iiwest+(jjwest-1)*npccx;
%                            ind = [indwestc indwestc+1];
%                            assert(max(ind)<=N,'index too large'); 
%                            Z(indf,ind) = 0.5; 
%          
%                        elseif(mod(i,2)==1 && mod(j,2)==0)
%                            %interpolate using south and north coarse grid
%                            %neighbors
%                            iisouth = ceil(i/2); jjsouth=ceil(j/2);
%                            indsouthc = iisouth+(jjsouth-1)*npccx;
%                            ind = [indsouthc indsouthc+npccx];                       
%                            assert(max(ind)<=N,'index too large'); 
%                            Z(indf,ind) = 0.5;
%                            
%                        elseif(mod(i,2)==0 && mod(j,2)==0)
%                            %interpolate using SW,SE,NW,NE coarse grid
%                            %neighbors
%                            iiSW = ceil(i/2);  jjSW = ceil(j/2);
%                            iiSE = iiSW+1;      jjSE = jjSW;
%                            iiNW = iiSW;        jjNW = jjSW+1;
%                            iiNE = iiSE;        jjNE = jjSW+1;
%                            ii = [iiSW iiSE iiNW iiNE];
%                            jj = [jjSW jjSE jjNW jjNE];
% 
%                            ind = ii+(jj-1)*npccx;
%                            assert(max(ind)<=N,'index is too large');
%                            Z(indf,ind) = 0.25;
%                        end
%                    end
%                end