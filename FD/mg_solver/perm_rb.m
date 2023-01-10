function [ P ] = perm_rb(s)
%% PERM_RB Constructs a permutation matrix for the red black ordering
%        on a square domain
%
%   Input:
%       s:    number of points in domain (must be square)
%   Output       
%       P:    Red-black permutation matrix
%             P*(1:n)'= (R1 R2 R3 R...B1 B2 B3...)';

%   Author: Luis Garcia Ramos,          
%           Institut fur Mathematik, TU Berlin
%           Version 1.0, Jun 2016
%
%
%%
assert(sqrt(s)==round(sqrt(s)),'input must be a square number'); 
%n: number of points in one dimension
n = sqrt(s);
I = speye(s);
%We construct the permutation matrix for the red-black ordering
if(mod(n,2)==1);
     k=round((s-1)/2);       
     P = sparse(s,s);
     %red variables
     P(:,1:2:(s+1)) = I(:,1:1:(k+1));
     %black variables
     P(:,2:2:(s-1)) = I(:,k+2:1:s);     
else   
    k=round(n/2);
    for i=1:k        
        v1 = ((2*(i-1)*n+1):2:2*(i-1)*n+n)';
        v2 = ((2*i-1)*n+2:2:(2*i-1)*n+n)';
        %red variables
        perm(n*(i-1)+1:n*i,1) = [v1;v2];       
        %black variables
        perm(2*k^2+n*(i-1)+1:2*k^2+n*i,1) = [v1+1;v2-1];       
    end    
    P = I(perm,:);      
end

end

