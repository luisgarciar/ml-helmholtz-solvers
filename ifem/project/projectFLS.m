%% Project: Fast Least Square Methods
%
% In this project we will implement fast least square methods by random
% sampling of a matrix transformed by randomized Fourier type transforms.
% We refer details of implementation to the first paper and theoretical
% background (original description) of algorithms to the rest:
%
% # Avron, H., Maymounkov, P., Toledo, S.: Blendenpik: Supercharging
% LAPACK?s least-squares solver. SIAM J. Scientific Comput. 32(3),
% 1217 - 1236, (2010).
% # Rokhlin, V., Tygert, M.: A fast randomized algorithm for overdetermined
% linear least-squares regression. Proc. Natl. Acad. Sci. USA 105(36),
% 13212 - 13217, (2008).
% # Drineas, P., Mahoney, M. W., Muthukrishnan, S., Sarlós, T. Faster least
% squares approximation. Numerische Mathematik, 117(2), 219 - 249, (2011).

%% Algorithm
%
% BLENDENPIK's Algorithm in p.1224 of [1].
%
% * Use deiscrete cosine transform (DCT) or discrete Hartley transform
% (DHT). 
% * set gamma = 4.
% * set the tolerence of LSQR as 10^{-6}

%% Experiments
%
% Test three types of matrices: X, Y and Z listed in p.1227 of [1].
% 
% * Robust to the condition number. Reproduce Fig. 5.4.
% * Cost of the different phases. Reproduce Fig. 5.7.
% * Convergence rate. Reproduce Fig. 5.6.