function [z] = pnebi1(x)
%FUNCTION y = pnebi1(x). Pnebi stands for (double of) product of negative 
% exponential (with absolute value of argument) and BesselI function i.e.,
% 2.0 * $\exp{-||x||)} \besseli_1(x)
% where Bessel function of order 1 is taken into consideration
% Note that Bessel functions are even
%

ax = abs(x);
  if (ax < 3.75) 
    y = (x/3.75)^2;
    z = ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934+y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
    z = 2.0 * exp(-ax) * z;
  else
    y=3.75/ax;
    z=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1-y*0.420059e-2));
    z=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2+y*(0.163801e-2+y*(-0.1031555e-1+y*z))));
    z = 2.0 * z / sqrt(ax);
  end
  if (x < 0.0)
     z = -z;
  end

end %file function


% #include <math.h>
% float bessi1(float x)
% Returns the modified Bessel function I1 (x) for any real x.
% {
% float ax,ans;
% double y;
% Accumulate polynomials in double precision.
% if ((ax=fabs(x)) < 3.75) {
% Polynomial fit.
% y=x/3.75;
% y*=y;
% ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
% +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
% } else {
% y=3.75/ax;
% ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
% -y*0.420059e-2));
% ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
% +y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
% ans *= (exp(ax)/sqrt(ax));
% }
% return x < 0.0 ? -ans : ans;
% }