classdef NonlinearSolver < handle
	% Abstract framework for the 

	properties (SetAccess = immutable)
		NonLinearArgs
	end
	
	properties (SetAccess = protected)
		LinearSolver
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
		end
		
		function obj = SetLinearSolver(obj,linsolve)
			obj.LinearSolver = linsolve;
		end
	end

	methods(Abstract)
		
		[xf, optout, stats] = solve(nonlinsystem, f, mass_scale, jac_scale, x0, optin);
	end
end

