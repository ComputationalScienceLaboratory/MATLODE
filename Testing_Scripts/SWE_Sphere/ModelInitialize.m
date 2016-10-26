%  Solves the SWE on the sphere problem from 
%  B. Neta, F.X. Giraldo and I.M. Navon
%  Analysis of the Turkel-Zwas Scheme for the Two-Dimensional 
%  Shallow Water Equations in Spherical Coordinates.

% Load the initialconditions. The initialconditions can also be built.
%See the BuildInitialConditions.m
load initialconditions
y0 = initialconditions;

% Suggested integration information.(Required)
Tspan = [0 1200]; % Corresponds to 1 Hr.
dt_integration = 50;


% RHS and Jacobian function calls, in MATLAB standard form.(Required)
rhsFun = @SWE_SphereTZ;
jacFun = @SWE_SphereJacSparse;

[tout, yout] = ode45(rhsFun, Tspan, y0);