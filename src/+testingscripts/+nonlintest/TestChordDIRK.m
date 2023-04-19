clear
format long e
close all

integrator = matlode.rk.dirk.ESDIRK_4_3_6;

options.ErrNorm = matlode.errnorm.InfNorm(1e-6, 1e-6);
options.StepSizeController = matlode.stepsizecontroller.StandardController;

problem = otp.kpr.presets.Canonical;

sol = integrator.integrate(problem.RHS, problem.TimeSpan, problem.Y0, options);

sol_matlab = problem.solve('Solver', @ode45, 'RelTol', 1e-6, 'AbsTol', 1e-6);

true_sol = problem.solve;

norm(sol.y(:,end) - true_sol.y(:,end))
norm(sol_matlab.y(:,end) - true_sol.y(:,end))
norm(sol.y(:,end) - sol_matlab.y(:,end))
