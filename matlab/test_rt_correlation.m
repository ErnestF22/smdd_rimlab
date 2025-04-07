function test_rt_correlation


N = 10;
points = rand(2, 2*N); % points vector
points(1,1:N) = 5.0;
points(2,1:N) = 1:N;
points(1,N+1:2*N) = 7.0;
points(2,N+1:2*N) = 1:N;
sigma = 0.5;

T.R = eye(2);
T.t = zeros(2,1);

%% compute correlation matrix

corr = rt_correlation(points, sigma, 1.0, T);
disp("corr")
disp(corr)

Ts = -2:0.1:5;

y = zeros(1, length(Ts));
progr_id = 1;
for t = Ts
    T.t(2) = t;
    corr = rt_correlation(points, sigma, 1.0, T);
    y(progr_id) = corr;
    disp("corr transformed")
    disp(corr)
    progr_id = progr_id + 1;
end

figure(1)
plot(Ts, y);

%%

theta = 0:0.01:2*pi;
% theta = pi;
rho = -12:0.1:12;
[rt, Theta, Rho] = rt_gmm_eval(points,sigma, 1.0, theta, rho);
figure(2)
contour(Theta, Rho, rt, 30)
% plot(rho, rt)

T.t(2) = 2.5;
pointsTransf = T.R * points + repmat(T.t, 1, 2*N);
figure(3)
[rt, Theta, Rho] = rt_gmm_eval(pointsTransf,sigma, 1.0, theta, rho);
contour(Theta, Rho, rt, 30)
% plot(rho, rt)

%% Correlation w.r.t. theta with T = identity
T.R = eye(2);
T.t = zeros(2,1);
T.t(1) = 0.6;
nf = 10;
nt = 1000;
[corr_coeffs] = rt_corr_fourier_coeffs(nf, points, sigma, T);
theta = 2*pi/nt * (0:nt-1)';
CosSin = zeros(nt, 2*nf+2);
for k=0:nf
  CosSin(:,2*k+1)   = cos(2 * k * theta);
  CosSin(:,2*k+2) = sin(2 * k * theta);
end
corr_theta = CosSin * corr_coeffs;
figure(4)
title('Fourier in T=eye(3)')
plot(theta, corr_theta);


end %file function