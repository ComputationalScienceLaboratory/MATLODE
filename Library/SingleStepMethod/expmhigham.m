function [E,normA] = expmhigham(A)
% EXPMHIGHAM - Matrix exponential via Pade approximation.
%   Part of Algorithm 2.3 in "The scaling and squaring method for
%   the matrix exponential revisited" by N. Higham, SIAM
%   J. Matrix Anal. Appl. 2005 (26) 1179 - 1193

% This is the first version need to add the case of using lower
% degree Pade approximations if the norm is less than the crucial
% theta_m quantity.
%
% JN: Changed to provide normA

% Scale A by power of 2 so that its norm is < 5.3.
  normA = norm(A, 1);
  tm = 5.4;
  [t s] = log2(normA / tm);
  s = max(0, s - (t == 0.5)); % adjust s if normA/theta(end) is a power of 2.
  A = A / 2^s;

  % Coefficients of degree 13 Pade approximation
  b = [64764752532480000, 32382376266240000, 7771770303897600, ...
       1187353796428800, 129060195264000, 10559470521600, 670442572800, ...
       33522128640, 1323241920,40840800, 960960, 16380, 182, 1];

  % Pade approximation for exp(A)
  A2 = A*A;
  A4 = A2*A2;
  A6 = A2*A4;

  U = A*(A6*(b(14)*A6+b(12)*A4+b(10)*A2)+b(8)*A6+b(6)*A4+b(4)*A2+b(2)* ...
         eye(size(A)));
  V = A6*(b(13)*A6+b(11)*A4+b(9)*A2)+b(7)*A6+b(5)*A4+b(3)*A2+b(1)* ...
         eye(size(A));

  E = (-U+V)\(U+V);

  % Undo scaling by repeated squaring
  for k=1:s, E = E*E; end
end
