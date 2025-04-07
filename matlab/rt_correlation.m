function corr = rt_correlation(points, sigma, weights, T)
%FUNCTION RT_CORRELATION() computes Radon Transform correlation for
% point set points, transform T (struct with R, t members) and variance sigma
%

N = size(points,2);

if isscalar(weights)
    weights = weights * ones(1,N);
end

%% compute correlation matrix

k = 1 / (2 * sqrt(pi) * abs(sigma)); %so simple due to isotropy

corr = 0.0;
for ii = 1:N
    for jj = 1:N
        M_ij = norm(T.R * points(:,jj) + T.t - points(:,ii));
%         phi_ij = (T.R * p(:,jj) + T.t - p(:,ii));
        pn0 = pnebi0(M_ij / (8 * sigma^2));
        corr = corr + pn0;
    end
end

corr = k * corr;


end %file function
