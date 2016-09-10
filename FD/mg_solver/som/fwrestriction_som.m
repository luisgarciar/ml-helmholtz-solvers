function R = fwrestriction_som(npf,dim,bc)
%% FWRESTRICTION_SOM Constructs the matrix corresponding to the full weight restriction
%   operator from a fine grid with npf interior points.
%  
%   Use:    R = fwrestriction(npf,dim,bc)  
%
%   Input: 
%       npf:    number of interior points in 1-D fine grid (must be odd)
%       dim:    dimension (1 or 2)
%       
%   Output:
%       R:      restriction matrix of size npc x npf
%
%   Author: Luis Garcia Ramos, 
%           Institut fur Mathematik, TU Berlin
%  
%  Version 1.0, Jun 2016
%  Works on 1-D, 2-D Dirichlet boundary conditions
%               
%%
switch dim
    case 1
        %npc = round(npf/2)-1; %length of coarse grid vectors
        y = zeros(npf,1); y(1:3,1) = [1;2;1];
        R = gallery('circul',y');
        R = 0.25*sparse(R(1:2:(npf-2),:));
           
    case 2
        switch bc
            case 'dir'
                %npc = round(npf/2)-1;  
                y   = zeros(npf,1); y(1:3,1) = [1;2;1];
                R   = gallery('circul',y)';
                R   = 0.25*sparse(R(:,1:2:(npf-2))); %1D operator%
                R   = kron(R,R)';  %2D operator  
                
            case 'som'
                assert(mod(npf,2)==1,'number of interior points must be even')
                npc = round((npf+1)/2)-1;
                npff  = npf+2; % total number of points with endpoints
                npcc  = npc+2;
                
   
                R     = sparse(npcc^2,npff^2);
                size(R)

                %The restriction matrix is filled by rows 
                %(change this later!)
                %indc: row index (coarse grid)
                %indf: column index (fine grid)

                %(0,0)- South-West Corner
                indc=1; indf=1; 
                R(indc,indf)=4; R(indc,indf+1)=4;
                R(indc,indf+npff)=4; R(indc,indf+npff+1)=4;
                
                %(1,0)-  South-East Corner
                indc=npcc; indf=npff;
                R(indc,indf)=4;  R(indc,indf-1)=4;
                R(indc,indf+npff)=4; R(indc,indf+npff-1)=4;
                
                %(0,1)- North-West Corner
                indc=npcc*(npcc-1)+1; indf=npff*(npff-1)+1;
                R(indc,indf)=4;   R(indc,indf-npf)=4;
                R(indc,indf+1)=4; R(indc,indf-npf+1)=4;
                
                %(1,1)- North East Corner
                indc=npcc^2; indf=npff^2;
                R(indc,indf)=4; R(indc,indf-1)=4;
                R(indc,indf-npff)=4; R(indc,indf-npff-1)=4;
                
                
                %South boundary y=0
                for indc=2:(npcc-1)
                    indf=2*indc-1;
                    R(indc,indf)=4; R(indc,indf+npff)=4;
                    R(indc,indf+1)=2; R(indc,indf-1)=2;
                    R(indc,indf+npff+1)=2; R(indc,indf+npff-1)=2; 
                end
                
                %East boundary x=1
                for i=2:(npcc-1)
                    indc=i*npcc; indf=(2*i-1)*npff;
                    R(indc,indf)=4; R(indc,indf-1)=4;
                    R(indc,indf-1+npff)=2; R(indc,indf-1-npff)=2;
                    R(indc,indf-npff)=2; R(indc,indf+npff)=2;
                end
                
                
                %North boundary y=1
                for i=2:(npcc-1)
                   indc=(npcc-1)*npcc+i; indf=(npff-1)*npff+(2*i-1);
                   R(indc,indf)=4; R(indc,indf-npf)=4; 
                   R(indc,indf+1)=2; R(indc,indf-1)=2; 
                   R(indc,indf-npff+1)=2; R(indc,indf-npff-1)=2;
                end
                
                %West boundary x=0
                for i=2:(npcc-1)
                   indc=npcc*(i-1)+1; indf=npff*2*(i-1)+1;
                   R(indc,indf)=4; R(indc,indf+1)=4;
                   R(indc,indf+npff)=2;  R(indc,indf-npff)=2;
                   R(indc,indf+1+npff)=2;  R(indc,indf+1-npff)=2;        
                end
                             
                %Interior points
                for i=2:(npcc-1)
                    for j=2:(npcc-1)
                        indc=i+npcc*(j-1);
                        ii=2*i-1; jj=2*j-1;
                        indf=ii+npff*(jj-1);
                        
                        R(indc,indf)=4;
                        R(indc,indf+npf)=2; R(indc,indf-npf)=2;
                        R(indc,indf+1)=2; R(indc,indf-1)=2;
                        R(indc,indf-1-npf)=1; R(indc,indf-1+npf)=1;
                        R(indc,indf+1-npf)=1; R(indc,indf+1+npf)=1;
                        
                    end
                end
                R=R/16;
                
        end      
      
end

end

