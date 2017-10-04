function ydot = vanDerPol_Function(~, y)

mu = 10.0;

% ydot = [y(2); mu*(1-y(1)^2)*y(2)-y(1)];
ydot = [y(2); -y(1) - mu*y(2)*(y(1)^2-1) ]; % <---using simple

return;

