clear
format long e
close all

%% Test Problem

A11 = [[-1, 2];...
	   [2, -1]];
A12 = [[2,0];...
	   [-1,-1]];
A21 = [[2,2];...
	   [4,4]];
A22 = [[-2,1];...
	   [1,-2]];

A = [[A11, A12];
	 [A21, A22]];

Mass = eye(4);
Mass(3,3) = 0;
Mass(4,4) = 0;

f = @(t,y) A*y;

tspan = [0,1];

Q = -(A22 \ A21);
S = A11 + A12 * Q;

x0 = [1;1];

ydiffExact = @(t) expmv(S, x0, t);
yalgExact = @(t) Q*expmv(S, x0, t);

yfullExact = @(t) [ydiffExact(t); yalgExact(t)];

z0 = yalgExact(tspan(1));

y0 = [x0;z0];


%% Integrator Parameters

integrator = matlode.rosenbrock.RODAS3;

options.ErrNorm = matlode.errnorm.InfNorm(1e-8, 1e-8);
options.StepSizeController = matlode.stepsizecontroller.StandardController;
% options.LinearSolver = matlode.linearsolver.MatrixLinearSolver(@mldivide, {'Jacobain', test_problem.jac});
options.Jacobian = A;
options.Mass = Mass;


%% Integration

uni_steps = linspace(0,1,160);

% sol = integrator.integrateFixed(f, uni_steps, y0, options);

sol = integrator.integrate(f, tspan, y0, options);

norm(sol.y(:,end) - yfullExact(tspan(end)))