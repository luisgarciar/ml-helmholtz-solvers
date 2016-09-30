%test of linear interpolation
 clear all;
 npc = 31;
 npcx  = npc; npcy = npc;
 npccx = npcx+2; npccy= npcy+2;
 npffx = 2*npccx-1; npffy = 2*npccy-1; 
 intdiv= @(x,y) idivide(int32(x),int32(y),'ceil');
 M = npffx*npffy; N = npccx*npccy;
 Z = sparse(M,N);
 
 npffx
 npffy
 npccx
 npccy
 
 size(Z)
 
 for j=1:npffx
     for i=1:npffy
         indf = i+(j-1)*npffx;          
        if (mod(i,2)==1 && mod(j,2)==1)
           %if fine grid point coincides with coarse grid point
           %transfer the function value unchanged
            %sequential index of point
           ii = intdiv(i,2); jj=intdiv(j,2); 
           indc = ii+(jj-1)*npccx;
           Z(indf,indc)=1;

        elseif (mod(i,2)==0 && mod(j,2)==1)
            %interpolate using east and west neighbors from
            %coarse grid
            iiwest = intdiv(i,2); jj=intdiv(j,2);
            %index of west neighbor in coarse grid
            indwestc = iiwest+(jj-1)*npccx; 
            ind = [indwestc indwestc+1];
            assert(max(ind)<=N,'index too large'); 
            Z(indf,ind) = 0.5; 

        elseif(mod(i,2)==1 && mod(j,2)==0)
               %interpolate using south and north coarse grid
               %neighbors
               iisouth   = intdiv(i,2); jjsouth=intdiv(j,2);
               indsouthc = iisouth+(jjsouth-1)*npccx;
               ind = [indsouthc indsouthc+npccx];
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
            ind = ii+(jj-1)*npccx;
            assert(max(ind)<=N,'index is too large');
            Z(indf,ind) = 0.25;
        end
     end
 end
 
 
 