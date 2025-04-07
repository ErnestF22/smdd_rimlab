function test_rt_gmm_eval

N = 10;
points = rand(2, 2*N); % points vector
points(1,1:N) = 5.0;
points(2,1:N) = 1:N;
points(1,N+1:2*N) = 7.0;
points(2,N+1:2*N) = 1:N;
sigma = 0.5;

rt_gmm_eval(points, sigma); %RADON



end %file function
