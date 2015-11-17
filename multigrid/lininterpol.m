function v_fine = lininterpol( v_coarse )
%INTERPOL Interpolates a coarse grid function to a fine grid in 1-D
%         using linear interpolation.
%   INPUT: 
%       v_coarse:    coarse grid function 
%       
%   OUTPUT:
%       v_fine:      fine grid function

%lengths of fine and coarse grid vectors
nc  = length(v_coarse);
nf  = 2*(nc+1)-1;

%interpolation
v_fine           = zeros(nf,1);
v_fine(2:2:nf)   = v_coarse;
v_fine(1:2:nf-2) = .5*(v_coarse);
v_fine(3:2:nf)   = v_fine(3:2:nf) +.5*(v_coarse); 

end

