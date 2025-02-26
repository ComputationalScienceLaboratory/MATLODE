classdef Chord < matlode.nonlinearsolver.NonlinearSolver
	%CHORD Classic Chord method for solving nonlinear systems. Does not
	%calculate the jacobian more than once
	properties(Access = protected)
		q0
		t0_pre
	end
	
	methods
		function obj = Chord(linsolve, args)
			arguments
				linsolve(1,1) matlode.linearsolver.LinearSolver = matlode.linearsolver.MatrixLinearSolver();
				args(1,:) cell = {};
			end
            
            obj = obj@matlode.nonlinearsolver.NonlinearSolver(linsolve, args{:});
		end


		function [out_opts, stats] = preprocess(obj, f, t0, y0, optin, stats)
			% Preprocess to compute M(t_0) y_0
			[stats] = obj.LinearSolver.computeMass(f, t0, y0, stats);
			if isempty(f.Mass) || ~isa(f.Mass, 'function_handle')
				obj.q0 = 0;
			else
				obj.q0 = (obj.LinearSolver.mass * y0);
			end
			obj.t0_pre = t0;
			out_opts = [];
		end
		
		function [xn, out_opts, stats] = solve(obj, f, t, dt, y0, x0, fn0, sys_const, mass_scale, jac_scale,  ~, stats)
			%Solve the nonlinear equation 0 = -M(t)z + const + a * f(t, y_0 + z) 
			xn = x0;
			i = 0;
			y1 = y0 + x0;
			intial_fun_cond = uint32(~isempty(fn0));
			if intial_fun_cond == 0
				fn0 = f.F(t, y1);
			end

			stats = obj.LinearSolver.preprocess(f, t, y1, true, -mass_scale, jac_scale, stats);

			while i < obj.MaxIterations
				mx = obj.LinearSolver.mass * x0;
				b =  (-mass_scale) .* (mx) + sys_const + (jac_scale) .* fn0;
				[w_i, stats] = obj.LinearSolver.solve(-b, stats);

				xn = w_i + x0;
				
				i = i + 1;

				if norm(w_i) < obj.Tolerance
					break;
				end
				x0 = xn;
				y1 = y0 + x0;
				fn0 = f.F(t, y1);
			end
			stats.nNonLinIterations = stats.nNonLinIterations + i;
			stats.nFevals = stats.nFevals + i - intial_fun_cond;

			out_opts.convergenceFailure = i >= obj.MaxIterations;
		end
	end

	methods (Access=private)
		
	end
end

