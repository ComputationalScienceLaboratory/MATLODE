function H = lorenz_Hess_vec( ~,~,u,v )

    H(1,1) = 0;
    H(2,1) = -v(3)*u(1) - v(1)*u(3);
    H(3,1) = v(2)*u(1) + v(1)*u(2);

return;

