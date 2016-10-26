%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           	    COMPUTATIONAL SCIENCE LABORATORY                   %%%
%%%                       MOFSET: Initialization                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Website: http://csl.cs.vt.edu/                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Suggested integration information.
tspan = [0 50];


m1 = 1;
m2 = 20;
k1 = 20;
k2 = 1;
k3 = 1;

x0 = [0; 0; 1; 0];
% RHS and Jacobian function calls, in MATLAB standard form.
rhsFun = @(t, x) springs_Function(t, x, m1, m2, k1, k2, k3);
%jacFun = @(t, U, Up) mosfet_Jacobian(t, U, Uop, U0, UT, fTilde);