classdef MatrixFreeLinearSolver < matlode.linearsolver.LinearSolver

    methods
        function obj = MatrixFreeLinearSolver(solver, solverArgs)
			arguments
				solver(1,1) function_handle = @gmres;
				solverArgs(1,:) cell = {};
			end
       
            obj = obj@matlode.linearsolver.LinearSolver(solver, solverArgs{:});
        end
        
        function [stats] = preprocess(obj, f, t, y, reeval, mass_scale, jac_scale, stats)
			if reeval
				if ~isempty(f.JacobianVectorProduct)
                	obj.jac = @(v) f.JacobianVectorProduct(t, y, v);
					stats.nJacobianEvals = stats.nJacobianEvals + 1;
				else
					if f.JLinearOperatorType ~= matlode.LinearOperatorType.Identity && f.JLinearOperatorType ~= matlode.LinearOperatorType.Constant
                		obj.jac = @(v) f.Jacobian(t, y) * v;
						stats.nJacobianEvals = stats.nJacobianEvals + 1;
					elseif isempty(obj.jac)
                		obj.jac = @(v) f.Jacobian * v;
					end
				end


				if ~isempty(f.MVectorProduct)
					obj.mass = @(v) f.MVectorProduct(t, y, v);
					stats.nMassEvals = stats.nMassEvals + 1;
				elseif f.MLinearOperatorType ~= matlode.LinearOperatorType.Empty
					if f.MLinearOperatorType ~= matlode.LinearOperatorType.Identity && f.MLinearOperatorType ~= matlode.LinearOperatorType.Constant
                		obj.mass = @(v) f.Mass(t, y) * v;
						stats.nMassEvals = stats.nMassEvals + 1;
					else
						if isempty(obj.mass)
							obj.mass = @(v) f.Mass * v;
						end
					end
				else
					if isempty(obj.mass)
						obj.mass = @(v) speye(size(y,1)) * v;
					end
				end
			end

            obj.system = @(v) mass_scale * obj.mass(v) + jac_scale * obj.jac(v);
        end
        
        function [sol, stats] = solve(obj, x, stats)
            [sol, flag, relres, iter] = obj.Solver(obj.system, x, obj.SolverArgs{:});
            
            stats.nLinearSolves = stats.nLinearSolves + 1;
        end
    end
end

