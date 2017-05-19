function H = vanDerPol_Hess_vec( ~,y,u,v )

    mu = 10.0;
    
    H(1,1) = 0;
    H(2,1) = -2*mu*(u(1)*(y(2)*v(1)+y(1)*v(2))+y(1)*v(1)*u(2));


return;

