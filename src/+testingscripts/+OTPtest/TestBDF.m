clear
format long e
close all

integrator = matlode.lmm.BDF();

problem = otp.allencahn.presets.Canonical;

options.ErrNorm = matlode.errnorm.InfNorm(1e-6, 1e-6);
options.StepSizeController = matlode.stepsizecontroller.StandardController;

base = 1e3;
tspan = (base .^ (linspace(0, 1, 200 + 1)) - 1) / (base - 1) * (problem.TimeSpan(end) - problem.TimeSpan(1)) + problem.TimeSpan(1);

sol = integrator.integrateFixed(problem.RHS, tspan, problem.Y0, options);
sol.stats

sol_matlab = problem.solve('RelTol', 1e-8, 'AbsTol', 1e-8);

norm(sol.y(end) - sol_matlab.y(end))
