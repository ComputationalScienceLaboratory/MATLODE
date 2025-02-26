classdef AdaptiveChord < matlode.nonlinearsolver.NonlinearSolver
	%An adaptive nonlinear solver utilizing Chord. If the problem doesn't
	%converge fast enough it will fail and ask for an easier problem

	properties(Access = protected)
		t0_pre
	end
	
	methods
		function obj = AdaptiveChord(linsolve, args)
			arguments
				linsolve(1,1) matlode.linearsolver.LinearSolver = matlode.linearsolver.MatrixLinearSolver();
				args(1,:) cell = {};
			end
            
            obj = obj@matlode.nonlinearsolver.NonlinearSolver(linsolve, args{:});
		end


		function [out_opts, stats] = preprocess(obj, f, t0, y0, mass_scale, jac_scale, optin, stats)
			% Preprocess the system

			stats = obj.LinearSolver.preprocess(f, t0, y0, true, -mass_scale, jac_scale, stats);

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
			out_opts.convergenceFailure = false;

			%Cannot assume that the system has updated mass to use for RHS
			%in the time depedent case
			[stats] = obj.LinearSolver.computeMass(f, t, y1, stats);

			nuRate = 2;
			deltaZ = 0;
			deltaZold = 0;


			while i < obj.MaxIterations
				mx = obj.LinearSolver.mass * x0;
				b =  (-mass_scale) .* (mx) + sys_const + (jac_scale) .* fn0;
				[w_i, stats] = obj.LinearSolver.solve(-b, stats);

				xn = w_i + x0;

				deltaZ = norm(w_i);

				if i > 0
					thetaRate = deltaZ ./ deltaZold;
					if thetaRate < 0.99
						nuRate = thetaRate ./ (1.0 - thetaRate);
						errPred = deltaZ .* thetaRate.^(obj.MaxIterations - i) ./  (1 - thetaRate);

						if errPred >= obj.Tolerance
							out_opts.convergenceFailure = true;
							break;
						end
					else
						out_opts.convergenceFailure = true;
                    	break;
					end
				end

                i = i + 1;
                
				if nuRate * deltaZ < obj.Tolerance
                    break;
				end

				deltaZold = deltaZ;
				x0 = xn;
				y1 = y0 + x0;
				fn0 = f.F(t, y1);
			end
			stats.nNonLinIterations = stats.nNonLinIterations + i;
			stats.nFevals = stats.nFevals + i - intial_fun_cond;

			out_opts.convergenceFailure = (out_opts.convergenceFailure || i >= obj.MaxIterations);
		end
	end

	methods (Access=private)
		
	end
end

