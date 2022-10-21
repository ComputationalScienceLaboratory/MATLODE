classdef MatrixFreeLinearSolver < matlode.linearsolver.LinearSolver

    methods
        function obj = MatrixFreeLinearSolver(solver, solverArgs)
			arguments
				solver(1,1) function_handle = @gmres;
				solverArgs(1,:) cell = {};
			end
       
            obj = obj@matlode.linearsolver.LinearSolver(solver, solverArgs{:});
        end
        
        function [stats] = preprocess(obj, t, y, reeval, mass_scale, jac_scale, stats)
			if reeval
                obj.mass = @(v) f.MVectorProduct(t, y, v);
				stats.nMassEvals = stats.nMassEvals + 1;

                obj.jac = @(v) f.JVectorProduct(t, y, v);
				stats.nJacobianEvals = stats.nJacobianEvals + 1;
			end

            obj.system = @(v) mass_scale * obj.mass(v) + jac_scale * obj.jac(v);
        end
        
        function [sol, stats] = solve(obj, x, stats)
            sol = obj.Solver(obj.system, x, obj.SolverArgs{:});
            
            stats.nLinearSolve = stats.nLinearSolve + 1;
        end
    end
end

