function jac = springs_Jacobian(~, ~, m1, m2, k1, k2, k3)

jac = [
    0 1 0 0;
    -(k1 + k2)/m1 0 k2/m1 0;
    0 0 0 1;
    k2/m2 0 -(k2+k3)/m2 0
];

end