classdef Chord < matlode.nonlinearsolver.NonlinearSolver
	%CHORD Classic Chord method for solving nonlinear systems. Does not
	%calculate the jacobian more than once
	
	properties
		Property1
	end
	
	methods
		function obj = Chord(linsolve, args)
			arguments
				linsolve(1,1) matlode.linearsolver.LinearSolver = matlode.linearsolver.MatrixLinearSolver();
				args(1,:) cell = {};
			end
            
            obj = obj@matlode.nonlinearsolver.NonlinearSolver(linsolve, args{:});
		end
		
		function [xn, xnf, optout, stats] = solve(obj, f, t, x0, sys_const, mass_scale, jac_scale,  ~, stats)
			xn = x0;
			i = 0;
			stats = obj.LinearSolver.preprocess(f, t, xn, true, mass_scale, jac_scale, stats);
			if isa(obj.LinearSolver.mass,'function_handle')
				mx = obj.LinearSolver.mass(xn);
			else
				mx =  obj.LinearSolver.mass * xn;
			end

			while i < obj.MaxIterations
				if i ~= 0
					[stats] = obj.LinearSolver.computeMass(f, t, xn, stats);
					if isa(obj.LinearSolver.mass,'function_handle')
						mx = obj.LinearSolver.mass(xn);
					else
						mx =  obj.LinearSolver.mass * xn;
					end
				end


				xnf = f.F(t,xn);
				b = -mx .* (mass_scale) + sys_const + jac_scale * xnf;
				[w_i, stats] = obj.LinearSolver.solve(b, stats);

				xtn = w_i + xn;
				i = i + 1;

				if norm(w_i) < obj.Tolerance
					break;
				end
				xn = xtn;
			end
			stats.nNonLinIterations = stats.nNonLinIterations + i;
			stats.nFevals = stats.nFevals + i;


			if i == obj.MaxIterations
				%TODO Add Flags
				optout = false;
			else
				optout = true;
			end
		
		end
	end

	methods (Access=private)
		
	end
end

