clear
format long e
close all

%% Test Problem

lambda = -1;
test_problem = testingscripts.testproblems.ODE.EulerProblem(lambda);

t_f = 1;

%% Integrator Parameters

integrator = matlode.rosenbrock.LinBackEuler;

options.ErrNorm = matlode.errnorm.InfNorm(1e-6, 1e-6);
options.StepSizeController = matlode.stepsizecontroller.StandardController;
options.LinearSolver = matlode.linearsolver.MatrixLinearSolver([], 'Jacobain', test_problem.jac);

%% Integration

sol = integrator.integrate(@(t,y) test_problem.f(t,y), [test_problem.t_0, t_f], test_problem.y_0, options);

norm(sol.y(end) - test_problem.y_exact(t_f))