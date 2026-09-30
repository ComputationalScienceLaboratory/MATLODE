clear
format long e
close all

integrator = matlode.lmm.BDF();

% options.ErrNorm = matlode.errnorm.InfNorm(1e-12, 1e-12);
% options.StepSizeController = matlode.stepsizecontroller.StandardController;

lambda = -1;
test_problem = testingscripts.testproblems.ODE.EulerProblem(lambda);

options.Jacobian = test_problem.jac;

t_f = 1;

sol = integrator.integrateFixed(@(t,y) test_problem.f(t,y), linspace(test_problem.t_0, t_f, 101), test_problem.y_0, options);

norm(sol.y(end) - test_problem.y_exact(t_f))