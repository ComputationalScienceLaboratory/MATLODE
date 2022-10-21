classdef MatrixLinearSolver < matlode.linearsolver.LinearSolver
    
    
    methods
        function obj = MatrixLinearSolver(solver, solverArgs)
			arguments
				solver(1,1) function_handle = @mldivide;
				solverArgs(1,:) cell = {};
			end
            
            obj = obj@matlode.linearsolver.LinearSolver(solver, solverArgs{:});
        end
        
        function [stats] = preprocess(obj, f, t, y, reeval, mass_scale, jac_scale, stats)
            if reeval
				if f.MLinearOperatorType ~= matlode.LinearOperatorType.Identity && f.MLinearOperatorType ~= matlode.LinearOperatorType.Constant
                	obj.mass = f.Mass(t, y);
					stats.nMassEvals = stats.nMassEvals + 1;
				elseif isempty(obj.mass)
					obj.mass = f.Mass;
				end

				if f.JLinearOperatorType ~= matlode.LinearOperatorType.Identity && f.JLinearOperatorType ~= matlode.LinearOperatorType.Constant
                	obj.jac = f.Jacobian(t, y);
					stats.nJacobianEvals = stats.nJacobianEvals + 1;
				elseif isempty(obj.jac)
                	obj.jac = f.Jacobian;
				end
            end

            obj.system = mass_scale * obj.mass + jac_scale * obj.jac;
        end
        
        function [sol, stats] = solve(obj, x, stats)
            sol = obj.Solver(obj.system, x, obj.SolverArgs{:});
            
            stats.nLinearSolves = stats.nLinearSolves + 1;
        end
    end
end

