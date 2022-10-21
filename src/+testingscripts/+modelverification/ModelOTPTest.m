
clear all
close all
format long e

% To test the interactions between MATLODE, OTP, and standard matlab
% options


problem = otp.lorenz63.presets.Canonical;

model = matlode.Model(problem.RHS);

mass = 1;

model2 = matlode.Model(problem.RHS, 'Mass', mass);

otpsets = odeset('Mass', mass, 'RelTol', 1e-6, 'AbsTol', 1e-6);
dfdt = @(t,y) 0;

model3 = matlode.Model(problem.RHS, otpsets, 'PartialDerivativeTime', dfdt);

