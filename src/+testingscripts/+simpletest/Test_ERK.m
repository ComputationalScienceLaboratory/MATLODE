clear
format long e
close all

integrator = matlode.rk.erk.DormandPrince();

options.ErrNorm = matlode.errnorm.InfNorm(1e-12, 1e-12);
options.StepSizeController = matlode.stepsizecontroller.StandardController;

lambda = -1;
test_problem = testingscripts.testproblems.ODE.EulerProblem(lambda);

t_f = 1;

sol = integrator.integrate(@(t,y) test_problem.f(t,y), [test_problem.t_0, t_f], test_problem.y_0, options);

norm(sol.y(end) - test_problem.y_exact(t_f))