function [coeffs] = rt_corr_fourier_coeffs(n,points, sigma, T)
% function [coeffs] = rt_corr_fourier_coeffs(n,points, sigma, T)
% Computes the Fourier coefficients of the correlation of pencil of
% parallel lines, i.e.
%  C(T) = \int RT[f](theta, rho) * RT[T[f]](theta, rho) drho
  np = size(points,2);
  if size(points,1) ~= 2
    disp(['points must have dimension np * 2, np is ' num2str(np)]);
    error('failure!')
  end
  
  pointsTransf = T.R * points + repmat(T.t,1,np);
  coeffs = zeros(2*n+2,1);
  norm_const = 1 / (2 * sqrt(pi) * abs(sigma)); %so simple due to isotropy
  for ii = 1:np
    for jj = 1:np
      dpoint = pointsTransf(:,jj) - points(:,ii);
      M_ij = norm(dpoint);
      phi_ij = atan2(dpoint(2),dpoint(1));
      pnebiM = norm_const .* pnebi(0:n, M_ij^2 / (8 * sigma^2));
      cosPhi = cos(2 .* (0:n) .* phi_ij) .* (-1).^(0:n);
      sinPhi = sin(2 .* (0:n) .* phi_ij) .* (-1).^(0:n);
      %size_prod = size((pnebiM .* cosPhi)')
      %size_coeff_odd = size(coeffs(1:2:(2*n+2)))
      %size_coeff_even = size(coeffs(2:2:(2*n+2)))
      coeffs(1:2:(2*n+2)) = coeffs(1:2:(2*n+2)) + (pnebiM .* cosPhi)';
      coeffs(2:2:(2*n+2)) = coeffs(2:2:(2*n+2)) + (pnebiM .* sinPhi)';
    end
  end
end