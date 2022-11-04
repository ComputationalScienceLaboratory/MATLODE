clear
format long e
close all

integrator = matlode.rosenbrock.RODAS4;

problem = otp.allencahn.presets.Canonical;

options.ErrNorm = matlode.errnorm.InfNorm(1e-8, 1e-8);
options.StepSizeController = matlode.stepsizecontroller.StandardController;
options.LinearSolver = matlode.linearsolver.MatrixFreeLinearSolver;

tic

sol = integrator.integrate(problem.RHS, problem.TimeSpan, problem.Y0, options);

toc

sol_matlab = problem.solve('Solver', @ode23t ,'RelTol', 1e-8, 'AbsTol', 1e-8);

norm(sol.y(end) - sol_matlab.y(end))
