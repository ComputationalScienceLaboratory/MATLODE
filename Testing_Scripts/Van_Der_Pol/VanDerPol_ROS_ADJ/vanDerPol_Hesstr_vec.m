function tmp = vanDerPol_Hesstr_vec( ~, y, u, k )
    t = 0;
    mu = 10;
    
    tmp(1,1) = -2*mu*u(2)*(y(2)*k(1)+y(1)*k(2));
    tmp(2,1) = -2*mu*u(2)*y(1)*k(1);

return;

