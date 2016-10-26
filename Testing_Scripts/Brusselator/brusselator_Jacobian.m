function jac = brusselator_Jacobian(~, ~, a, b)
%VDPOLJAC van der Pol equation Jacobian.
    
jac(1,1) = b - 1;
jac(1,2) = a;
jac(2,1) = -b;
jac(2,2) = -a;

end