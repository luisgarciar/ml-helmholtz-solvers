function v_coarse = fwrestriction( v_fine )
%FWRESTRICTION Restricts a fine grid function to a coarse grid in 1-D
%using full-weight restriction.
%   INPUT: 
%       v_fine:      fine grid function
%       
%   OUTPUT:
%       v_coarse:    coarse grid function 


%lengths of fine and coarse grid vectors
nf = length(v_fine);
nc = (.5*(nf+1))-1;

%restriction
v_coarse = .5*v_fine(2:2:nf);
v_coarse = v_coarse + .25*(v_fine(1:2:nf-2))+.25*(v_fine(3:2:nf));

end

