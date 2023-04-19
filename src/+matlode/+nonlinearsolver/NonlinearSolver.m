classdef NonlinearSolver < handle
	% Abstract framework for the 

	properties (SetAccess = immutable)
		NonLinearArgs
	end
	
	properties (SetAccess = protected)
		LinearSolver
		MaxIterations
		Tolerance
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
		
		[xn, xnf, optout, stats] = solve(obj,f, t, x0, sys_const, mass_scale, jac_scale,  optin, stats);
	end
end

