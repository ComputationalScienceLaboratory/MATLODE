clear
format long e
close all

integrator = matlode.rk.dirk.SDIRK_2_1_2;

options.ErrNorm = matlode.errnorm.InfNorm(1e-8, 1e-8);
options.StepSizeController = matlode.stepsizecontroller.StandardController;
options.NonLinearSolver = matlode.nonlinearsolver.Newton;

problem = otp.robertson.presets.Canonical;

sol = integrator.integrate(problem.RHS, problem.TimeSpan, problem.Y0, options);

sol_matlab = problem.solve('Solver', @ode15s, 'RelTol', 1e-8, 'AbsTol', 1e-8);

true_sol = problem.solve;

norm(sol.y(:,end) - true_sol.y(:,end))
norm(sol_matlab.y(:,end) - true_sol.y(:,end))
norm(sol.y(:,end) - sol_matlab.y(:,end))
