function Z = lininterpol(npc,dim,bc)
%% LININTERPOL Constructs the matrix corresponding to the (bi)linear 
%interpolation operator from a coarse grid with npc interior points.
%  
%   Use:    Z = lininterpol_som(npc,dim,bc)  
%
%   Input: 
%       npc:  number of 1D interior points in coarse grid
%       dim:  dimension (1 or 2)
%       bc:  'dir' for homogeneous Dirichlet boundary conditions
%            'som' for first order Sommerfeld boundary conditions 
%       
%   Output:
%       Z:  interpolation matrix, size npf x npc     (1D Dirichlet)
%                                      npf^2 x npc^2 (2D Dirichlet)
%                                     (npf+2)^2 x (npc+2)^2 (2D Sommerfeld)
%
%  Author: Luis Garcia Ramos, 
%          Institut fur Mathematik, TU Berlin
%          Version 2.0, Jun 2016
%
%  Dirichlet and Sommerfeld boundary conditions
%
%  References: 
%  A Multigrid Tutorial, Briggs, Henson, McCormick, 2000, Chap 3
%  Computational Science and Engineering, Strang, 2007, Chap 7 
%  Multigrid, Trottenberg, Oosterlee, Schiller, 2001, Chap 2
%
%
%%

switch dim
    case 1
        npf = 2*(npc)+1; %length of fine grid vectors
        z   = zeros(npf,1); z(1:3,1) = [1;2;1];
        Z   = gallery('circul',z)';
        Z   = 0.5*sparse(Z(:,1:2:(npf-2)));
        
    case 2
        npf = 2*(npc)+1; %number of interior points in fine grid (1D)
        switch bc
            case 'dir'                
                z   = zeros(npf,1); z(1:3,1) = [1;2;1];
                Z   = gallery('circul',z)';
                Z   = 0.5*sparse(Z(:,1:2:(npf-2))); %1D operator
                Z   = kron(Z,Z);  %2D operator 
                
            case 'som'               
               %modify this lines later to add option for
               %diff. number of points on x and y
               npcix  = npc; npciy = npc;
               npcx = npcix+2; npcy= npciy+2;
               npfx = 2*npcx-1; npfy = 2*npcy-1;
               intdiv = @(x,y) idivide(int32(x),int32(y),'ceil');
                             
               M = npfx*npfy; N = npcx*npcy;
               Z=sparse(M,N);

               %We fill the matrix by rows
               for j=1:npfy
                   for i=1:npfx
                       %(i,j): lexicographic index of point (hx*(i-1),hy*(j-1))
                       %in fine grid
                       indf = i+(j-1)*npfx; %sequential index of point
                       if (mod(i,2)==1 && mod(j,2)==1)
                       %if fine grid point coincides with coarse grid point
                       %transfer the function value unchanged
                           ii = intdiv(i,2); jj = intdiv(j,2);
                           ii = intdiv(i,2); jj=intdiv(j,2); 
                           indc = ii+(jj-1)*npcx;
                           assert(indc<=N,'index too large')
                           Z(indf,indc) = 1;
                           
                       elseif (mod(i,2)==0 && mod(j,2)==1)
                           %interpolate using east and west coarse grid
                           %neighbors 
                           iiwest = intdiv(i,2); jjwest=intdiv(j,2);
                           %index of west neighbor in coarse grid
                           indwestc = iiwest+(jjwest-1)*npcx;
                           ind = [indwestc indwestc+1];
                           assert(max(ind)<=N,'index too large'); 
                           Z(indf,ind) = 0.5; 
         
                       elseif(mod(i,2)==1 && mod(j,2)==0)
                           %interpolate using south and north coarse grid
                           %neighbors
                           iisouth = intdiv(i,2); jjsouth=intdiv(j,2);
                           indsouthc = iisouth+(jjsouth-1)*npcx;
                           ind = [indsouthc indsouthc+npcx];                       
                           assert(max(ind)<=N,'index too large'); 
                           Z(indf,ind) = 0.5;
                           
                       elseif(mod(i,2)==0 && mod(j,2)==0)
                           %interpolate using SW,SE,NW,NE coarse grid
                           %neighbors
                           iiSW = intdiv(i,2); jjSW = intdiv(j,2);
                           iiSE = iiSW+1;      jjSE = jjSW;
                           iiNW = iiSW;        jjNW = jjSW+1;
                           iiNE = iiSE;        jjNE = jjSW+1;
                           ii = [iiSW iiSE iiNW iiNE];
                           jj = [jjSW jjSE jjNW jjNE];

                           ind = ii+(jj-1)*npcx;
                           assert(max(ind)<=N,'index is too large');
                           Z(indf,ind) = 0.25;
                       end
                   end
               end
        end
end
        
end

