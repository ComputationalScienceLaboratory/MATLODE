classdef NonlinearSolver < handle
	% Abstract framework for the 

	properties (SetAccess = immutable)
		NonLinearArgs
	end
	
	properties (SetAccess = protected, GetAccess = public)
		LinearSolver
		MaxIterations
		Tolerance
		AbsTol
		RelTol
		RatioTol
	end
	
	methods
		function obj = NonlinearSolver(linsolve, args, maxiter, tolerance)
			arguments
				linsolve(1,1) matlode.linearsolver.LinearSolver = {};
				args(1,:) cell = {};
				maxiter(1,1) int64 = 100
				tolerance(1,1) double = 1e-8
			end
            
            obj.NonLinearArgs = args;
			obj.LinearSolver = linsolve;
			obj.MaxIterations = maxiter;
			obj.Tolerance = tolerance;
		end
		
		function obj = SetLinearSolver(obj,linsolve)
			obj.LinearSolver = linsolve;
		end
	end

	methods(Abstract)

		[out_opts, stats] = preprocess(obj, f, t0, y0, mass_scale, jac_scale, optin, stats);
		
		[xn, xnf, out_opts, stats] = solve(obj, f, t, dt, y0, x0, fn0, sys_const, mass_scale, jac_scale,  optin, stats);
	end
end

