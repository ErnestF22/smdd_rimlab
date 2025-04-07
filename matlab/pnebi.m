function [y] = pnebi(n,x)
%FUNCTION y = pnebi(n,x). Pnebi stands for (double of) product of negative 
% exponential (with absolute value of argument) and BesselI function i.e.,
% 2.0 * $\exp{-||x||)} \besseli_n(x)
% where Bessel function of order k is taken into consideration
% Note that Bessel functions are even
%
   BIG_NUM = 1e+10;
   SMALL_NUM = 1e-10;
   x = abs(x);    % PNEBI is an even function
   N = max(n);    % n can be a vector of indices, e.g. n=1:6, so N is the max, e.g. 6
   % PNEBI(k,0) = 0 for k > 0, PNEBI(0,0) = 2
   Y = zeros(N+1);
   Y(1) = 2.0;
   if x < SMALL_NUM
     y = Y(n+1);
     return;
   end
   M = 2 * (N + ceil(sqrt(40.0 * N)));
   Y = zeros(N+1);
   factor = 2.0 / x;
   seqPrev = 0.0;  % bip
   seqCurr = 1.0;  % bi
   seqNext = 0.0;  % bim
   for k=M:-1:0
     %disp('k:');
     %disp(k)
     seqNext = seqPrev + factor * k * seqCurr;
     seqPrev = seqCurr;
     seqCurr = seqNext;
     % It only stores in Y the values of PNEBI in interval 0:N
     if k <= N 
       Y(k+1) = seqPrev;
     end
     % To avoid overflow!
     if seqCurr > BIG_NUM
       seqPrev = seqPrev * SMALL_NUM;
       seqCurr = seqCurr * SMALL_NUM;
       Y = Y * SMALL_NUM;
       %disp('Rescale to avoid overflow')
     end
   end
   scaleFactor = pnebi0(x) / Y(1);
   Y = scaleFactor * Y;
   y = Y(n+1);
end
