function y_p = arenstorfOrbit_Function( t, y )
%%% Allocate Memory %%%
y_p = zeros(4,1);

%%% INPUT PARAMETERS %%%
% y(1) = q1 = x1
% y(2) = r3 = x2
% y(3) = q2 = x1_p
% y(4) = r4 = x2_p
mu = 0.012277471;
mu_hat = 1 - mu;

%%% OUTPUT ARGUMENTS %%%
% f(1) = u1_pp
% f(2) = u2_pp

%D1 = ( ( y(1)+mu )*( y(1)+mu ) + y(2)*y(2) )^(1.5);
%D2 = ( ( y(1)-mu_hat )*( y(1)-mu_hat ) + y(2)*y(2) )^(1.5);

%f(1) = y(1) + 2*y(4) - mu_hat*( y(1)+mu )/D1 - mu*( y(1)-mu_hat )/D2;
%f(2) = y(2) - 2*y(3) - mu_hat*y(2)/D1 - mu*y(2)/D2;
%f(3) = D1;
%f(4) = D2;

%%% NOT USING SIMPLE() %%%
% y_p(1) = y(3);
% y_p(2) = y(4);
% y_p(3) = y(1) + 2*y(4) - mu_hat*( y(1)+mu )/( ( ( y(1)+mu )^2 + y(2)^2 )^(3/2) ) ...
%     - mu*( y(1)-mu_hat )/( ( ( y(1)-mu_hat )^2 + y(2)^2 )^(3/2) );
% y_p(4) = y(2) - 2*y(3) - mu_hat*y(2)/( ( ( y(1)+mu)^2 + y(2)^2 )^(3/2) ) ...
%     - mu*y(2)/( ( ( y(1)-mu_hat )^2 + y(2)^2 )^(3/2) );

%%% USING SIMPLE %%%
y_p(1,1) = y(3);
y_p(2,1) = y(4);
y_p(3,1) = y(1) + 2*y(4) - (mu*(mu + y(1) - 1))/((mu + y(1) - 1)^2 + y(2)^2)^(3/2) + ((mu + y(1))*(mu - 1))/((mu + y(1))^2 + y(2)^2)^(3/2);
y_p(4,1) = ((mu - 1)/((mu + y(1))^2 + y(2)^2)^(3/2) - mu/((mu + y(1) - 1)^2 + y(2)^2)^(3/2) + 1)*y(2) - 2*y(3);

return;

