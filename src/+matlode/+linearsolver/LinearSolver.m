classdef (Abstract) LinearSolver < handle
    
    properties (SetAccess = immutable, GetAccess = protected)
        Solver
		SolverArgs
    end

    properties (SetAccess = protected)
		mass
		jac
        system
	end
    
    methods
        function obj = LinearSolver(solver, solverArgs)
			arguments
				solver(1,1) function_handle;
				solverArgs(1,:) cell = {};
			end
            
            obj.Solver = solver;
			obj.SolverArgs = solverArgs;
        end
    end
    
    methods (Abstract)
        % mass_scale*mass + jac_scale*jac
        [stats] = preprocess(obj, f, t, y, reeval, mass_scale, jac_scale, stats);

		[stats] = computeMass(obj, f, t, y, stats);
        
        [sol, stats] = solve(obj, x, stats);
        
        
    end
end

