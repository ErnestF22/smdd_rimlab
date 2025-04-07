function [y] = pnebi0(x)
%FUNCTION y = pnebi0(x). Pnebi stands for (double of) product of negative 
% exponential (with absolute value of argument) and BesselI function i.e.,
% 2.0 * $\exp{-||x||)} \besseli_0(x)
% where Bessel function of order 0 is taken into consideration
% Note that Bessel functions are even
%

  x = abs(x);
  t = x / 3.75;
  if (t < 1.0) 
    t2 = t*t;
    y = 1.0 + t2 * (3.5156229 + t2 * (3.0899424 + t2 * (1.2067492 + t2 * (0.2659732 + t2 * (0.360768e-1 + t2 * 0.45813e-2)))));
    y = 2.0 * exp(-x) * y;
  else
    tinv = 1 / t;
    y = (0.39894228 + tinv * (0.1328592e-1 + tinv * (0.225319e-2 + tinv * (-0.157565e-2 + tinv * ...
        (0.916281e-2 + tinv * (-0.2057706e-1 + tinv * (0.2635537e-1 + tinv * (-0.1647633e-1 + tinv * 0.392377e-2))))))));
    y = 2.0 * y / sqrt(x);
  end

 end %file function