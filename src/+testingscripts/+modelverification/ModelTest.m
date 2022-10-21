
clear all
close all
format long e

%Test should ignore RelTol AbsTol without error
%Test will copy model 1 into model 2
%Test should modify model and set operator types

f = @(t,y) y;
jac = @(t,y) 1;
mass = 1;
dfdt = @(t,y) 0;

model = matlode.Model(f, 'Jacobian', jac, 'Mass', mass, 'PartialDerivativeTime', dfdt, 'RelTol', 1e-6, 'AbsTol', 1e-6);

f2 = @(t,y) 2*y;

model2 = matlode.Model(model);

jac2 = @(t,y) t;
mass2 = 2;

otpsets = odeset('Jacobian', jac2, 'Mass', mass2, 'RelTol', 1e-6, 'AbsTol', 1e-6);

model3 = matlode.Model(model, otpsets);