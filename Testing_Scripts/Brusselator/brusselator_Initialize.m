%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           	    COMPUTATIONAL SCIENCE LABORATORY                   %%%
%%%                     Van Der Pol: Initialization                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Website: http://csl.cs.vt.edu/                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Suggested integration information.
tspan = [0 10];

a = 1;
b = 4.5;

% linear representation of initial condition.
x0 = [2; 1];

% RHS and Jacobian function calls, in MATLAB standard form.
rhsFun = @(t, y) brusselator_Function(t, y, a, b);
jacFun = @(t, y) brusselator_Jacobian(t, y, a, b);