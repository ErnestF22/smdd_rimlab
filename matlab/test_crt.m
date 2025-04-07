function test_crt
 
  N = 10;
  points = rand(2, 2*N); % points vector
  points(1,1:N) = 5.0;  
  points(2,1:N) = 1:N;
  points(1,N+1:2*N) = 7.0;  
  points(2,N+1:2*N) = 1:N;
  sigma = 0.5;
  weight = 1.0;
  T.R = eye(2);
  T.t = zeros(2,1);
  T.t(1) = 5.0;
  pointsTransf = T.R * points + repmat(T.t, 1, 2*N);

  [mus_out,sigmas_out,weights_out] = crt_computation(points, sigma, weight);
  disp("mus_out size:")
  disp(size(mus_out))
  disp("sigmas_out size:")
  disp(size(sigmas_out))
  disp("weights_out size:")
  disp(size(weights_out))

  theta = 0:0.01:2*pi;
  rho = -12:0.1:12;

  [rt1, Theta, Rho] = rt_gmm_eval(points,sigma,1.0, theta, rho);
  figure(1)
  contour(Theta, Rho, rt1, 30)
  title("RT")

  [rt2, Theta, Rho] = rt_gmm_eval(pointsTransf,sigma,1.0, theta, rho);
  figure(2)
  contour(Theta, Rho, rt2, 30)
  title("RT Transformed")

  [crt1, Theta, Rho] = rt_gmm_eval(mus_out,sigmas_out,weights_out, theta, rho);
  figure(3)
  %contour(Theta, Rho, crt1, 30)
  surfc(Theta, Rho, crt1)
  title("CRT")

  [mus_out,sigmas_out,weights_out] = crt_computation(pointsTransf, sigma, weight);  
  [crt2, Theta, Rho] = rt_gmm_eval(mus_out,sigmas_out,weights_out, theta, rho);
  figure(4)
  %contour(Theta, Rho, crt2, 30)
  surfc(Theta, Rho, crt2)
  title("CRT Transformed")


  
end