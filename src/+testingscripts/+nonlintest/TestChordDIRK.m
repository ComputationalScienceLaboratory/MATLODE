clear
format long e
close all

integrator = matlode.rk.dirk.SDIRK_2_1_2;

options.ErrNorm = matlode.errnorm.StandardNorm(1e-7, 1e-7);
options.StepSizeController = matlode.stepsizecontroller.StandardController;
options.NonLinearSolver = matlode.nonlinearsolver.Chord;
options.NonLinearSolver.MaxIterations = 1000;

problem = otp.robertson.presets.Canonical;

tic
sol = integrator.integrate(problem.RHS, problem.TimeSpan, problem.Y0, options);
toc

tic
sol_matlab = problem.solve('Solver', @ode15s, 'RelTol', 1e-8, 'AbsTol', 1e-8);
toc

% beta = problem.Parameters.Beta;
% t = problem.TimeSpan(end);
% true_sol.y = [t .* sin(t) + (1 + beta * t) .* exp(-t); ...
%                 beta * exp(-t) + sin(t)];

true_sol = problem.solve( 'RelTol', 1e-10, 'AbsTol', 1e-10);
% true_sol.y = exp(-1);

norm(sol.y(:,end) - true_sol.y(:,end))
norm(sol_matlab.y(:,end) - true_sol.y(:,end))
norm(sol.y(:,end) - sol_matlab.y(:,end))

sol.stats
