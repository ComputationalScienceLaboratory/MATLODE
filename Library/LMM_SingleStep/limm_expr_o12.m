function [alpha, beta, gamma] = limm_expr_o12(k)
%
switch(k)
  case 1
    alpha = @(~) ...
      [(-1)];
    beta  = @(~) ...
      [1];
    gamma = @(~) ...
      [(-1),1];
  case 2
    alpha = @(w21) ...
      [0,(-1)];
    beta  = @(w21) ...
      [0,1];
    gamma = @(w21) ...
      [(1/2).*w21.^(-1),(-1)+(-1/2).*w21.^(-1),1];
  otherwise
    error('bad input k');
end
end
