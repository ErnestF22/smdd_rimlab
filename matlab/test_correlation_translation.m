function test_correlation_translation
  N = 10;
  points = rand(2, 2*N); % points vector
  points(1,1:N) = 5.0;    
  points(2,1:N) = 1:N;
  points(1,N+1:2*N) = 7.0;  
  points(2,N+1:2*N) = 1:N;
  sigma = 0.5;

  T.R = eye(2);
  T.t = zeros(2,1);

  transl_min = -5.0;
  transl_res = 0.1;
  transl_num = 100;
  corrRT = zeros(transl_num,transl_num);
  for ix=1:transl_num
    for iy=1:transl_num
      T.t(1) = transl_min + transl_res * ix;
      T.t(2) = transl_min + transl_res * iy;
      corrRT(ix,iy) = rt_correlation(points, sigma, 1.0, T);
    end
  end
  tval = transl_min + transl_res * (1:transl_num);
  [Tx,Ty] = meshgrid(tval,tval);
  figure(1)
  %contour(Tx,Ty,corrRT,50);
  surf(Tx,Ty,corrRT);
  title('Correlation of RT vs translation')
  xlabel('tx')
  ylabel('ty')
  zlabel('corr')

  [points_crt, sigma_crt, weights_crt] = crt_computation(points,sigma,1.0);
  corrCRT = zeros(transl_num,transl_num);
  for ix=1:transl_num
    for iy=1:transl_num
      T.t(1) = transl_min + transl_res * ix;
      T.t(2) = transl_min + transl_res * iy;
      corrCRT(ix,iy) = rt_correlation(points, sigma, weights_crt, T);
    end
  end
  tval = transl_min + transl_res * (1:transl_num);
  [Tx,Ty] = meshgrid(tval,tval);
  figure(2)
  %contour(Tx,Ty,corrRT,50);
  surf(Tx,Ty,corrCRT);
  title('Correlation of CRT vs translation')
  xlabel('tx')
  ylabel('ty')
  zlabel('corr')


end