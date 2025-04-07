function [low, mid, upp, level] = interval_pow2(v1, v2, res)
% FUNCTION [low, mid, upp] = intervalPow2(v1, v2, res)
% The function operates on discretized number i1 = round(v1/res) and i2 = round(v1/res) 
% where res is the resolution.
%
% OUTPUT:
% - low: maximum low = 2^l s.t. low <= i1 and low <= i2
% - upp: minimum upp = 2^u s.t. i1 <= upp and i2 <= upp
% - mid: the medium bewteen low and upp
% - level: minimum size of power 2 interval 2^level s.t. i1 and i2 are both
%   contained inside
  if (nargin < 3)
    res = 1;
  end
  i1 = round(v1 ./ res);
  i2 = round(v2 ./ res);
  if i1 == i2
    low = i1;
    mid = i1;
    upp = i1;
    level = 0;
    return;
  elseif (i1 > i2)
    tmp = i1;
    i1 = i2;
    i2 = tmp;
  end
  % e.g. in 8 bit representation, i1 00010110 and i2 00001101
  % level = floor(log2(00011011)) = 4
  % inverv_mask: 00000111
  % 
  level = floor(log2(bitxor(i1, i2))) + 1;
  lmax = ceil(log2(i2));
  interv_mask =  bitshift(1,level) - 1;
  negated_mask = bitxor(bitshift(1,lmax) - 1, interv_mask);
  low = bitand(i1, negated_mask);
  upp = bitor(low, interv_mask);
  mid = bitor(low, bitshift(1,level-1));
end