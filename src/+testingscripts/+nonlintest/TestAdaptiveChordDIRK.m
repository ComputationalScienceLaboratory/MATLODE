clear
format long e
close all

integrator = matlode.rk.dirk.SDIRK_2_1_2;

options.ErrNorm = matlode.errnorm.StandardNorm(1e-5, 1e-5);
options.StepSizeController = matlode.stepsizecontroller.StandardController;
options.NonLinearSolver = matlode.nonlinearsolver.Chord;

% problem = otp.robertson.presets.Canonical;



alpha = 1/2;

f = @(t, u) [[(alpha*(2-t) - 1)/((2-t)*(3-t)),1-t, (2-t)*alpha];[(alpha - 1)/((2-t)*(3-t)), -1, alpha-1]]  * u + [((3-t)/(2-t)); 2]*exp(t);
g = @(t,u) [(t+2)/(3-t), t^2 - 4,0] * u - (t^2 + t - 2)*exp(t);
f_full = @(t,u) [[(alpha*(2-t) - 1)/((2-t)*(3-t)),1-t, (2-t)*alpha];...
			  [(alpha - 1)/((2-t)*(3-t)), -1, alpha-1];...
			  [(t+2)/(3-t), t^2 - 4,0]] * u + [((3-t)/(2-t)); 2; - (t^2 + t - 2)]*exp(t);
g_u = @(t,u) [(t+2)/(3-t), t^2 - 4] * f(t,u) + [5/(t-3)^2, 2*t,0] * u - (2*t + 1)*exp(t) - (t^2 + t - 2)*exp(t);

jac = @(t,u) [[(alpha*(2-t) - 1)/((2-t)*(3-t)),1-t, (2-t)*alpha];...
			  [(alpha - 1)/((2-t)*(3-t)), -1, alpha-1];...
			  [(t+2)/(3-t), t^2 - 4,0]];

mass = [[1,0,0];...
		[0,1,0];...
		[0,0,0]];

jac_Q = @(t,u) [(t+2)/(3-t), t^2 - 4,0] * [(2-t)*alpha; alpha-1; 0];

t_0 = 0;
t_f = 1;
tspan = [t_0, t_f];

u_exact = @(t) [(3-t)*exp(t);exp(t);-exp(t) ./ (2-t)];

u_0 = u_exact(t_0);

dim_y = 2;


model = matlode.Model(f_full, 'Jacobian', jac, 'Mass', mass);

% tspan = linspace(t_0,t_f,2^9);


tic
% sol = integrator.integrate(problem.RHS, problem.TimeSpan, problem.Y0, options);
sol = integrator.integrate(model, tspan, u_0, options);
toc

% tic
% sol_matlab = problem.solve('Solver', @ode15s, 'RelTol', 1e-8, 'AbsTol', 1e-8);
% toc

% beta = problem.Parameters.Beta;
% t = problem.TimeSpan(end);
% true_sol.y = [t .* sin(t) + (1 + beta * t) .* exp(-t); ...
%                 beta * exp(-t) + sin(t)];
true_sol.y = u_exact(tspan(end));

% true_sol = problem.solve('Solver', @ode15s, 'RelTol', 1e-10, 'AbsTol', 1e-10);
% true_sol.y = exp(-1);

norm(sol.y(:,end) - true_sol.y(:,end))
% norm(sol_matlab.y(:,end) - true_sol.y(:,end))
% norm(sol.y(:,end) - sol_matlab.y(:,end))

sol.stats


