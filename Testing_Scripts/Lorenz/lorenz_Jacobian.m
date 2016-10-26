function [ jac ] = lorenz_Jacobian( ~, input_param, constant )
%LORENZ_JACOBIAN Summary of this function goes here
%   Detailed explanation goes here

    syms x y z
    x = input_param(1);
    y = input_param(2);
    z = input_param(3);
    
    if ( nargin < 3 )
        beta = 8/3;
        rho  = 28;
        sig  = 10;
    else
        beta = constant(1);
        rho  = constant(2);
        sig  = constant(3);
    end

    jac = [ -sig    sig     0; ...
            rho-z   -1      -x; ...
            y       x       -beta ];


return;

