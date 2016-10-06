function R = fwrestrictionfem(npf,dim)
%% FWRESTRICTIONFEM Constructs the matrix corresponding to the full weight restriction
%   operator from a fine grid with npf interior points.
%  
%   Use:    R = fwrestriction(npf,dim)  
%
%   Input: 
%       npf:    number of interior points in 1-D fine grid (must be even)
%       dim:    dimension (1 or 2)
%       
%   Output:
%       R:      restriction matrix of size npc x npf
%
%   Author: Luis Garcia Ramos, 
%           Institut fur Mathematik, TU Berlin
%  
%  Version 1.0, Jun 2016
%               
%%
switch dim
    case 1
        assert(mod(npf,2)==0,'number of gridpoints is not even');
        npc = round(npf/2); %length of coarse grid vectors
        y = zeros(npf,1); y(1:3,1) = [1;2;1];
        R = gallery('circul',y');
        R = 0.5*sparse(R(1:2:(npf),:));
        R(npc,:)= zeros(1,npf);R(npc,npf)=1; R(npc,npf-1)=0.5;
    case 2
%         %npc = round(npf/2)-1;  
%         y   = zeros(npf,1); y(1:3,1) = [1;2;1];
%         R   = gallery('circul',y)';
%         R   = 0.25*sparse(R(:,1:2:(npf-2))); %1D operator%
%         R   = kron(R,R)';  %2D operator   
%       
end

end

