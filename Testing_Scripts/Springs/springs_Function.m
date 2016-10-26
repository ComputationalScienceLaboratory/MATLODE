function xPrime = springs_Function(~, x, m1, m2, k1, k2, k3)

xPrime = [
    x(2);
    (-(k1 + k2) * x(1) + k2 * x(3)) / m1;
    x(4);
    (k2 * x(1) - (k2 + k3) * x(3)) / m2;
];

return;