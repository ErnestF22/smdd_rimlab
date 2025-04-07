function [mus_out,sigmas_out,weights_out] = crt_computation(mus_in,sigmas_in,weights_in)
  n_in = size(mus_in, 2);
  n_out = n_in * (n_in+1) / 2;
  
  if isscalar(sigmas_in)
    sigmas_in = sigmas_in * ones(1,n_in);
  end
  if isscalar(weights_in)
    weights_in = weights_in * ones(1,n_in);
  end

  progr_id = 0;
  size_mus = size(mus_in, 1);
  mus_out = zeros(size_mus, n_out);
  sigmas_out = zeros(1, n_out);
  weights_out = zeros(1, n_out);
  for ii=1:n_in
    for jj=ii:n_in
      sigma_sq_i = sigmas_in(ii)^2;
      sigma_sq_j = sigmas_in(jj)^2;
      sigma_sq_ij = (1 / sigma_sq_i + 1 / sigma_sq_j)^(-1);
      mu_ij = sigma_sq_ij .* (mus_in(:,ii) ./ sigma_sq_i +  mus_in(:,jj) ./ sigma_sq_j);
      sigma_sq_part = sigma_sq_i + sigma_sq_j;
      weight_part = exp(-0.5 * (mus_in(jj) - mus_in(jj))^2 / sigma_sq_part) / (sqrt(2*pi*sigma_sq_part));
      weight_ij = weights_in(ii) * weights_in(jj) * weight_part;
      if (jj > ii) 
          weight_ij = 2 * weight_ij;
      end
      progr_id = progr_id + 1;
      mus_out(:, progr_id) = mu_ij;
      sigmas_out(progr_id) = sqrt(sigma_sq_ij);
      weights_out(progr_id) = weight_ij;
    end
  end
  [weights_out, idx] = sort(weights_out);
  mus_out(:,idx) = mus_out(:,idx);
  sigmas_out = sigmas_out(idx);
end