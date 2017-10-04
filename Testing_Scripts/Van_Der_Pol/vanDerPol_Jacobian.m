function jac = vanDerPol_Jacobian(~, y)
%VDPOLJAC van der Pol equation Jacobian.

mu = 10.0;

% jac = [         0                    1
%        (-2*mu*y(1)*y(2)-1) (mu*(1-y(1)^2))];
    
jac(1,1) = 0;                   %
jac(1,2) = 1;                   %
                                % <--- using simple
jac(2,1) = -2*mu*y(1)*y(2)-1;   %
jac(2,2) = -mu*y(1)^2 + mu;     %

end