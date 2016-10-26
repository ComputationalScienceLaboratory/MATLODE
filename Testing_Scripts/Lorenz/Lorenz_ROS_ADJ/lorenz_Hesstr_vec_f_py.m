function tmp = lorenz_Hesstr_vec_f_py( ~, ~, u, k )

    tmp(1,1) = -u(3)*k(3);
    tmp(2,1) = u(2)*k(1);
    tmp(3,1) = -u(1)*(k(1)-k(2));

return;

