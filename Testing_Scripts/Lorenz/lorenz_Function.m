function ydot = lorenz_Function( ~, input_param, constant )

    x = input_param(1);
    y = input_param(2);
    z = input_param(3);

    if ( nargin < 3 )
        beta = 8/3;
        rho  = 28; % originallly 28 but unstable
        sig  = 10;
    else
        beta = constant(1);
        rho  = constant(2);
        sig  = constant(3);
    end

    ydot = [-sig*x+sig*y; rho*x-x*z-y; x*y-beta*z];

end

