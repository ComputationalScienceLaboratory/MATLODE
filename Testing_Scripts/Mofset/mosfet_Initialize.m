%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           	    COMPUTATIONAL SCIENCE LABORATORY                   %%%
%%%                       MOFSET: Initialization                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Website: http://csl.cs.vt.edu/                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Suggested integration information.
tspan = [0 50];

% The constant operating voltage
Uop = 5;

% The threshold voltage
UT = 1;

% The ground voltage
% Note this is not the vector representing the inital conditions
U0 = 0;

fTilde = @(UG, UD, US) max(UG - US - UT, 0).^2 - max(UG - UD - UT, 0).^2;

% linear representation of initial condition.
% All nodes start with 5 volts
Uinit = zeros(3,1) + 5;

% RHS and Jacobian function calls, in MATLAB standard form.
rhsFun = @(t, U, Up) mosfet_Function(t, U, Uop, U0, UT, fTilde);
jacFun = @(t, U, Up) mosfet_Jacobian(t, U, Uop, U0, UT, fTilde);