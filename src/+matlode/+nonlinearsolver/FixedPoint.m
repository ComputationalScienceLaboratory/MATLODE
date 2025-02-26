classdef FixedPoint < matlode.nonlinearsolver.NonlinearSolver
	%CHORD Classic Fixed Point iteration
	properties(Access = protected)
		t0_pre
	end
	
	methods
		function obj = FixedPoint(linsolve, args)
			arguments
				linsolve(1,1) matlode.linearsolver.LinearSolver = matlode.linearsolver.MatrixLinearSolver();
				args(1,:) cell = {};
			end
            
            obj = obj@matlode.nonlinearsolver.NonlinearSolver(linsolve, args{:});
		end


		function [out_opts, stats] = preprocess(obj, f, t0, y0, mass_scale, jac_scale, optin, stats)
			% Preprocess to compute M(t_0) y_0
			[stats] = obj.LinearSolver.computeMass(f, t0, y0, stats);
			%TODO: Fixed point can be formulated to support DAEs. Utilize
			%QR factorization plus some transformation to get it.
			if ~isempty(f.MassSingular)
				error('Fixed Point Iteration does not support DAEs.')
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

			stats = obj.LinearSolver.preprocess(f, t, y1, true, mass_scale, 0, stats);

			while i < obj.MaxIterations
				b = sys_const + (jac_scale) .* fn0;
				[xn, stats] = obj.LinearSolver.solve(b, stats);
				
				i = i + 1;

				if norm(xn - x0) < obj.Tolerance
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

