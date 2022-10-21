clear
format long e
close all

integrator = matlode.rk.erk.Vener87;

options.ErrNorm = matlode.errnorm.InfNorm(1e-8, 1e-8);
options.StepSizeController = matlode.stepsizecontroller.StandardController;

problem = otp.linear.presets.Canonical;

sol = integrator.integrate(problem.RHS, problem.TimeSpan, problem.Y0, options);

sol_matlab = problem.solve('RelTol', 1e-8, 'AbsTol', 1e-8);

true_sol = exp(-problem.TimeSpan(end));

norm(sol.y(end) - true_sol)
norm(sol_matlab.y(end) - true_sol)
norm(sol.y(end) - sol_matlab.y(end))
