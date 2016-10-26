function xdot = brusselator_Function(~, x, a, b)

temp = a * x(1)^2 * x(2);

xdot = [
    1 - (1 + b) * x(1) + temp;
    b * x(1) - temp
];

return;