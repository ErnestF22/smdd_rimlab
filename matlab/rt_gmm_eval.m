function [rt,Theta,Rho] = rt_gmm_eval(points, sigma, weight, theta, rho)
  % points is 2 * np
  % theta is nt * 1
  np = size(points,2);
  nt = size(theta,1);
  nr = size(rho, 1);

  if length(sigma) == 1
    sigmas = sigma * ones(1,np);
  else 
    sigmas = sigma;
  end
  if length(weight) == 1
    weights = weight * ones(1,np);
  else
    weights = weight;
  end
%   disp("size(sigmas)")
%   disp(size(sigmas))
%   disp("size(weights)")
%   disp(size(weights))
%   disp("np")
%   disp(np)
  assert(length(sigmas) == np & length(weights) == np)

  [Theta, Rho] = meshgrid(theta, rho);
  rt = zeros(nr,nt);
  for ii = 1:np
    aa = cos(Theta) * points(1,ii) + sin(Theta) * points(2,ii) - Rho;
    rt = rt + weights(ii) / sqrt(2 * pi * sigmas(ii)^2) * exp(-aa.^2 ./ (2 * sigmas(ii)));
  end
end