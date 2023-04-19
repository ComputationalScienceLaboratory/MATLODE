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
			%TODO: Add Case for array of t and y and matrix jac_scale for
			%FIRK
			if reeval


				if ~isempty(f.JacobianVectorProduct)
                	obj.jac = @(v) f.JacobianVectorProduct(t, y, v);
					stats.nJacobianEvals = stats.nJacobianEvals + 1;
				else
					if isa(f.Jacobian,'double')
						%TODO: Replace with Data Type
						if isempty(obj.jac)
							obj.jac = @(v) f.Jacobian * v;
						end
					else 
						obj.jac = @(v) f.Jacobian(t, y) * v;
						stats.nJacobianEvals = stats.nJacobianEvals + 1;
					end
				end

				stats = obj.computeMass(f, t, y, stats);
				
			end

            obj.system = @(v) mass_scale * obj.mass(v) + jac_scale * obj.jac(v);
		end

		function  [stats] = computeMass(obj, f, t, y, stats)
			if ~isempty(f.MVectorProduct)
				obj.mass = @(v) f.MVectorProduct(t, y, v);
				stats.nMassEvals = stats.nMassEvals + 1;
			else
				if isempty(f.Mass)
					if isempty(obj.mass)
						obj.mass = @(v) v;
					end
				elseif isa(f.Mass, 'double')
					%TODO: Replace with Data Type
					if isempty(obj.mass)
						obj.mass = @(v) f.Mass * v;
					end
				else
					obj.mass = @(v) f.Mass(t,y) * v;
					stats.nMassEvals = stats.nMassEvals + 1;
				end
			end
		end
        
        function [sol, stats] = solve(obj, x, stats)
            [sol, flag, relres, iter] = obj.Solver(obj.system, x, obj.SolverArgs{:});
            %TODO: Add Flag results
            stats.nLinearSolves = stats.nLinearSolves + 1;
        end
    end
end

