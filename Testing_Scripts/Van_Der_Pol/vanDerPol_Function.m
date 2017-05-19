function ydot = vanDerPol_Function(~,y, mu)

if nargin < 3 
    mu = 10.0;
end

% ydot = [y(2); mu*(1-y(1)^2)*y(2)-y(1)];
ydot = [y(2); -y(1) - mu*y(2)*(y(1)^2-1) ]; % <---using simple

return;

