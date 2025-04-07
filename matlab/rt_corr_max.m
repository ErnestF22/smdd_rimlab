function [corr] = rt_corr_max(points, sigma, T)
  nf = 20;
  nt = 200; 
  [corr_coeffs] = rt_corr_fourier_coeffs(nf, points, sigma, T);
  theta = pi/nt * (0:nt-1)';
  CosSin = zeros(nt, 2*nf+2);
  for k=0:nf
    CosSin(:,2*k+1) = cos(2 * k * theta);
    CosSin(:,2*k+2) = sin(2 * k * theta);
  end
  corr_theta = CosSin * corr_coeffs;
  corr = max(corr_theta);
end