function tmp = vanDerPol_Hesstr_vec_f_py( ~, y, u, k )

    tmp(1) = -2*u(2)*y(1)*y(2)*k(1) + u(2)*k(2)*(1-y(1)^2);

return;

