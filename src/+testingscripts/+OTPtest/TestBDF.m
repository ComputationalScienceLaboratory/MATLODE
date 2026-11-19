clear
format long e
close all

integrator = matlode.lmm.BDF();

problem = otp.lorenz63.presets.Canonical;

options.ErrNorm = matlode.errnorm.InfNorm(1e-6, 1e-6);
options.StepSizeController = matlode.stepsizecontroller.StandardController;
options.FullTrajectory = true;

base = 1e3;
tspan = (base .^ (linspace(0, 1, 200 + 1)) - 1) / (base - 1) * (problem.TimeSpan(end) - problem.TimeSpan(1)) + problem.TimeSpan(1);

% sol = integrator.integrateFixed(problem.RHS, tspan, problem.Y0, options);
sol = integrator.integrate(problem.RHS, problem.TimeSpan, problem.Y0, options);
sol.stats

sol_matlab = problem.solve('RelTol', 1e-8, 'AbsTol', 1e-8);

norm(sol.y(end) - sol_matlab.y(end)) / max(norm(sol.y(end)), norm(sol_matlab.y(end)))

f = figure;
t = tiledlayout('vertical')

nexttile
plot(sol.t, sol.y(1, :))

nexttile
plot(sol.t, sol.stats.order)

saveas(f, "order.svg")