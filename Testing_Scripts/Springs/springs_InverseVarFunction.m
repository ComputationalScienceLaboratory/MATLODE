function y = springs_InverseVarFunction(x)

w1 = x(1);
w2 = x(2);

C11 = x(3);
C13 = x(4);
C21 = x(5);
C23 = x(6);

a = C11 / C21;
b = C13 / C23;
denom = a * b * (w1^2 - w2^2);

k2 = 1;

m1 = k2 * (a - b) / denom;
m2 = k2 * a * b * (b - a) / denom;
k1 = k2 * (a * (1 - b) * w1^2 + b * (1 - a) * w2^2) / denom;
k3 = k2 * a * b * ((b - 1) * w1^2 - (1 - a) * w2^2) / denom;

x1 = C11 + C13;
x2 = C21 + C23;

y = [m1; m2; k1; k3; x1; x2] - [1; 20; 20; 1; 0; 1];

end