clear
format long e
close all

integrator = matlode.rk.erk.DormandPrince();

options.ErrNorm = matlode.errnorm.InfNorm(1e-12, 1e-12);
options.StepSizeController = matlode.stepsizecontroller.StandardController;

lambda = -1;
test_problem = testingscripts.testproblems.ODE.EulerProblem(lambda);

t_f = 1;

order_checker = matlode.util.OrderChecking(integrator);

stepamount = 2.^(5:8);
fig_num = 1;

[error, ~] = order_checker.getErrorWRTExactValue(@(t,y) test_problem.f(t,y), test_problem.t_0, t_f, test_problem.y_0, test_problem.y_exact(t_f), stepamount, options);

[f, polyError] = order_checker.plot_error_single(stepamount, error, fig_num)


sol = integrator.integrate(@(t,y) test_problem.f(t,y), [test_problem.t_0, t_f], test_problem.y_0, options);


norm(sol.y(end) - test_problem.y_exact(t_f))