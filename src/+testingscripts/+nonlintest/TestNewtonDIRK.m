clear
format long e
close all

integrator = matlode.rk.dirk.ESDIRK_2_1_3;

options.ErrNorm = matlode.errnorm.InfNorm(1e-7, 1e-7);
options.StepSizeController = matlode.stepsizecontroller.StandardController;
options.NonLinearSolver = matlode.nonlinearsolver.Newton;

problem = otp.robertson.presets.Canonical;

tic
sol = integrator.integrate(problem.RHS, problem.TimeSpan, problem.Y0, options);
toc

sol_matlab = problem.solve('Solver', @ode15s, 'RelTol', 1e-8, 'AbsTol', 1e-8);

true_sol = problem.solve;

norm(sol.y(:,end) - true_sol.y(:,end))
norm(sol_matlab.y(:,end) - true_sol.y(:,end))
norm(sol.y(:,end) - sol_matlab.y(:,end))
sol.stats
