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
			%TODO: Add Case for array of t and y and matrix jac_scale for
			%FIRK
            if reeval
				if isa(f.Jacobian,'double')
					%TODO: Replace with Data Type
					if  isempty(obj.jac)
						obj.jac = f.Jacobian;
					end
				else 
					obj.jac = f.Jacobian(t, y);
					stats.nJacobianEvals = stats.nJacobianEvals + 1;
				end

				stats = obj.computeMass(f, t, y, stats);
				
            end

            obj.system = mass_scale * obj.mass + jac_scale * obj.jac;
		end

		function [stats] = computeMass(obj, f, t, y, stats)
			if isempty(f.Mass)
				if isempty(obj.mass)
					if issparse(obj.jac)
						obj.mass = speye(length(y));
					else
						obj.mass = eye(length(y));
					end
				end
				
			elseif isa(f.Mass, 'double') 
				%TODO: Replace with Data Type
				if isempty(obj.mass)
					obj.mass = f.Mass;
				end
			else
				obj.mass = f.Mass(t,y);
				stats.nMassEvals = stats.nMassEvals + 1;
			end
		end
        
        function [sol, stats] = solve(obj, x, stats)
            sol = obj.Solver(obj.system, x, obj.SolverArgs{:});
            
            stats.nLinearSolves = stats.nLinearSolves + 1;
        end
    end
end

