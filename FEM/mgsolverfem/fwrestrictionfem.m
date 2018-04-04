function R = fwrestrictionfem(npf,dim,bc)
%% FWRESTRICTIONFEM Constructs the matrix corresponding to the full weight restriction
%   operator from a fine grid with npf interior points.
%  
%   Use:    R = fwrestriction(npf,dim,bc)  
%
%   Input: 
%       npf:    number of interior points in 1-D fine grid (must be odd)
%       dim:    dimension (1 or 2)
%       bc:     boundary conditions ('mix' or 'som')
%       
%   Output:
%       R:  restriction matrix
%           size(R) = ((npf+1)/2,npf+1)    if bc = 'mix'
%                   = ((npf-1)/2+2,npf+2) if bc = 'som'
%
%   Author: Luis Garcia Ramos, 
%           Institut fur Mathematik, TU Berlin
%  
%  Version 1.0, Jun 2016
%               
%%
switch dim
    case 1
        switch bc
            case 'mix'
                assert(mod(npf,2)==1,'number of interior gridpoints is not odd');
                npc = round((npf+1)/2); %length of coarse grid vectors
                y = sparse(npf+1,1); y(1:3,1) = [1;2;1];
                R = gallery('circul',y');
                %R = smcirc(y');
                R            = 0.5*sparse(R(1:2:(npf+1),:));
                R(npc,:)     = zeros(1,npf+1);
                R(npc,npf+1) = 1; R(npc,npf)=0.5;
            
            case 'som'
                assert(mod(npf,2)==1,'number of interior gridpoints is not odd');
                npc = round((npf-1)/2)+2; %length of coarse grid vectors
                y = zeros(npf+2,1); y(1:4,1) = [0;1;2;1];
                R1 = sparse(gallery('circul',y'));
                %R1 = smcirc(y');
                R1 = 0.5*sparse(R1(1:2:npf-1,:));
                R = sparse(npc,npf+2);
                R(2:npc-1,:) = R1;
                R(npc,npf+2)= 1; R(npc,npf+1)=0.5;
                R(1,1)=1; R(1,2)=0.5;
                
        end
                
                
    case 2
%         %npc = round(npf/2)-1;  
%         y   = zeros(npf,1); y(1:3,1) = [1;2;1];
%         R   = gallery('circul',y)';
%         R   = 0.25*sparse(R(:,1:2:(npf-2))); %1D operator%
%         R   = kron(R,R)';  %2D operator   
%       
end

end

