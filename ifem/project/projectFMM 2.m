%% Project: Fast Multipole Methods
%
% The purpose of this project is to write the tree code and fast multipole
% methods for the N-body summation problem.
%
% Reference: 
%
% * <http://math.uci.edu/~chenlong/226/FMMsimple.pdf Fast Multipole Methods>

%% Step 1: Direct Sum
%
% Generate two random vectors x, y with length N. Although the direct sum
% can be implemented in the double for loops, in MATLAB, it is better to
% generate the matrix first and then compute the matrix-vector product for
% another random vector q.

%% Step 2: Compute the weight
%
% To store the weight in different level, use |cell| structure. Use a for
% loop of i=1:N to compute the weight first and then try to remove this for
% loop to speed up your code. The loop over levels is small (only logN
% times) and thus can be kept.

%% Step 3: Evaluation
%
% First write code to find the interaction list. Then loop over each cell
% in a given level and compute the far field in the interaction list. In
% the fines level, add the near field by direct sum or matrix-vector
% product using a small matrix.
%

%% Step 4: Test
%
% # Chose N small and J = 1. Make sure the code works for one level (only
% four intervals) first by comparing the result using tree algorithm with
% the result in Step 1.
% # Test the performance for different N and plot the CPU time vs N for
% both direct method and tree code.

%% Step 5 (optional): Fast Multipole Methods
%
% Modify the tree code to fast multipole methods.
%
% # Compute the weight by restriction from the fine grid to coarse grid.
% # Implement the M2L: multipole expansion to local expansion
% # Change the evaluation of far field in the interaction list to the merge
% of coefficients b in the local expansion.
% # Translate the local expansion using the prolongation operator.
% # Evaluate in the finest level.
% # Plot the CPU time vs N to confirm the O(N) complexity.