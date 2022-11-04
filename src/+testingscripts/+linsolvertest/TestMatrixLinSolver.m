clear
format long e
close all

integrator = matlode.rosenbrock.ROS34PW2();

problem = otp.allencahn.presets.Canonical;

options.ErrNorm = matlode.errnorm.InfNorm(1e-6, 1e-6);
options.StepSizeController = matlode.stepsizecontroller.StandardController;
options.Jacobian = problem.RHS.Jacobian(problem.TimeSpan(1),problem.Y0);

tic

sol = integrator.integrate(problem.RHS, problem.TimeSpan, problem.Y0, options);

toc

sol_matlab = problem.solve('RelTol', 1e-12, 'AbsTol', 1e-12);

norm(sol.y(end) - sol_matlab.y(end))
