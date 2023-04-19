classdef Broyden < matlode.nonlinearsolver.NonlinearSolver
	%BROYDEN Classic Broyden Method. Works without a linear solver
	
	methods
		function obj = Broyden(args)
			arguments
				args(1,:) cell = {};
			end
            
            obj = obj@matlode.nonlinearsolver.NonlinearSolver([], args{:});
		end
		
		function [xf, optout, stats] = solve(obj, f, t, x0, sys_const, mass_scale, jac_scale,  ~, stats)
			xf = x0;
			i = 0;
			while i < obj.MaxIterations
				stats = obj.LinearSolver.preprocess(f, t, xf, true, mass_scale, jac_scale, stats);
				if isa(obj.LinearSolver.mass,'function_handle')
					mx = obj.LinearSolver.mass(xf);
				else
					mx =  obj.LinearSolver.mass * xf;
				end


				b = mx .* (-mass_scale) - sys_const - jac_scale * f.F(t,xf);
				[w_i, stats] = obj.LinearSolver.solve(b, stats);

				xtn = w_i + xf;
				i = i + 1;
				stats.nNonLinIterations = stats.nNonLinIterations + 1;

				if norm(w_i) < obj.Tolerance
					return;
				end
				xf = xtn;
			end

			if i == obj.MaxIterations
				%TODO Add Flags
				optout = false;
			end
		
		end
	end
end

