function test_rt_corr_max_manopt


N = 10;
points = rand(2, 2*N); % points vector
points(1,1:N) = 5.0;
points(2,1:N) = 1:N;
points(1,N+1:2*N) = 7.0;
points(2,N+1:2*N) = 1:N;
sigma = 0.5;

% Create the problem structure.
manifold = specialeuclideanfactory(2);
problem.M = manifold;

% Define the problem cost function and its Euclidean gradient.
problem.cost  = @(x) -rt_corr_max(points, sigma, x);
% problem.egrad = @(x) -2*A*x;      % notice the 'e' in 'egrad' for Euclidean
problem = manoptAD(problem); % automatic differentiation

% Numerically check gradient consistency (just once, optional).
checkgradient(problem);

options.maxiter = 50;

T_initguess.R = eye(2);
T_initguess.t = [0; -0.1];

disp("cost at initial guess")
disp(problem.cost(T_initguess))

T_eye.R = eye(2);
T_eye.t = zeros(2,1);
disp("cost at identity")
disp(problem.cost(T_eye))

% Solve.
[xOpt, xOptCost, info, options] = trustregions(problem, T_initguess, options);
disp("xOpt.R")
disp(xOpt.R)
disp("xOpt.t")
disp(xOpt.t)
disp("xOptCost")
disp(xOptCost)



end
