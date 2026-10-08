classdef BDF < matlode.Integrator
	properties (Constant)
		PartitionMethod = false;
		PartitionNum = 1;
		MultirateMethod = false;
	end

	properties
		NonLinearSolver
	end

	properties (SetAccess = protected)
		Pascal   % Pascal matrix needed for the predictor step. Same for all methods, should only vary in size based on MaxOrder.

		L        % The matrix of methods coefficients. Each row is a different order method, so there should be a number of rows equal to MaxOrder
		C        % Vector of error constants. Each entry is the error constant of a different order method
		MaxOrder % The maximum feasible order attainable by the method. Determines the size of the nordsieck vector.
	end

	methods
		function obj = BDF(datatype)
			arguments
				datatype(1,1) string = 'double';
			end

			obj = obj@matlode.Integrator(true, datatype);

			% In the future, these will be set by the individual method
			obj.MaxOrder = 5;
			obj.L = [...
				[1  , 1  , 0  , 0 , 0 , 0] / 1;...
				[2  , 3  , 1  , 0 , 0 , 0] / 3;...
				[6  , 11 , 6  , 1 , 0 , 0] / 11;...
				[24 , 50 , 35 , 10, 1 , 0] / 50;...
				[120, 274, 225, 85, 15, 1] / 274;...
				];
			obj.C = [-1/2, 2/9, -3/22, 12/125, -10/137];

			obj.Pascal = zeros(obj.MaxOrder + 1, obj.MaxOrder + 1);
			% Note that the pascal matrix constructed here is transposed from the text.
			% This is due to how we apply it to the nordsieck vector as a multivector
			for i = 1:obj.MaxOrder + 1
				obj.Pascal(i, 1) = 1;
				for j = 2:i
					obj.Pascal(i, j) = obj.Pascal(i - 1, j) + obj.Pascal(i - 1, j - 1);
				end
			end
		end
	end

	methods (Access = protected)
		function opts = matlodeSets(obj, p, varargin)

			%BDF Specific options
			p.addParameter('NonLinearSolver', matlode.nonlinearsolver.Chord());

			% TODO - add option for startup procedure

			opts = matlodeSets@matlode.Integrator(obj, p, varargin{:});

			if isempty(opts.NonLinearSolver)
				error('Please provide appropiate parameters and a non-linear solver')
			end
			obj.NonLinearSolver = opts.NonLinearSolver;
		end

		function [t, y, stats] = timeLoop(obj, f, tspan, y0, opts)

			numVars = length(y0);
			multiTspan = length(tspan) > 2;

			errNorm = opts.ErrNorm;
			stepController = opts.StepSizeController;

			tlen = 0;
			if multiTspan
				y = zeros(numVars, length(tspan));
				t = zeros(1, length(tspan));
			elseif opts.FullTrajectory
				y = zeros(numVars, opts.ChunkSize);
				t = zeros(1, opts.ChunkSize);
				tlen = length(t);
			else
				y = zeros(length(y0), 2);
				t = tspan;
			end

			t(:,1) = tspan(:,1);
			y(:, 1) = y0;

			tcur = tspan(1);
			tnext = tcur;
			dtbuf = 0; % Kahan summation buffer to remember small parts of dt which haven't yet been added to t
			tindex = 2;

			tspanlen = length(tspan);

			tdir = sign(tspan(end) - tspan(1));
			dtmax = min([abs(opts.MaxStep), abs(tspan(end) - tspan(1))]) * tdir;
			dtmin = opts.MinStep * tdir;

			%inital values
			ynext = y0;

			%Start stats
			stats = obj.intalizeStats;
			stats.nSteps = 1;

			%First step
			%allocates memory for first step
			% TODO - can get this from startup procedure
			order_next = 1;
			[dt0, f0, stats.nFevals] = stepController.startingStep(f, tspan, y0, order_next, errNorm, dtmin, dtmax);
			dtnext = dt0;
			dtprev = dt0;

			% nordsieck vector stored as a multivector, one column for each entry
			% TODO - get more values from startup procedure
			nordsieck_next = zeros(length(y0), obj.MaxOrder + 1);
			nordsieck_next(:, 1) = y0;
			nordsieck_next(:, 2) = dtnext * f0;

			prevAccept = true;

			% time loop
			while tindex <= tspanlen
				if(stats.nSteps > opts.MaxNumSteps)
					warning("OneStepIntegrator:MaxNumSteps", "Integrator reached maximum steps before finishing integration at t = %f.", tcur)
					break
				end

				ycur = ynext;
				tcur = tnext;
				dtprev = dtcur;
				dtcur = dtnext;
				nordsieck = nordsieck_next;
				order_current = order_next;

				% Accept Loop
				% Will keep looping until accepted step
				while true
					[ynext, nordsieck_next, delta, stats, out_opts] = obj.timeStep(f, tcur, ycur, dtcur, nordsieck, true, order_current, stats);
				
					if out_opts.failure == false
						[err, stats] = timeStepErr(obj, ycur, ynext, dtcur, errNorm, order_current, delta, stats, opts);

						% Find next Step
						% TODO: Update to take stats
						[prevAccept, dtnext] = stepController.newStepSize(prevAccept, tspan, dtcur, err, order_current);

						% Check the step size if we decrease order
						if order_current > 1
							% The nordsieck vector already stores an estimate of y_p, which is the residual for order p - 1 
							err = errNorm.errEstimate(ycur, ynext, obj.C(order_current - 1) * factorial(order_current) * nordsieck_next(:, order_current + 1));
							[~, dtnext_candidate] = stepController.newStepSize(prevAccept, tspan, dtcur, err, order_current - 1);

							if dtnext_candidate > dtnext
								dtnext = dtnext_candidate;
								order_next = order_current - 1;
							end
						end

						if order_current < obj.MaxOrder
							temp = factorial(current_order);
							temp2 = current_order * temp;
							% The p-1 derivative of y at last time step
							fnm1 = temp * nordsieck(:, order_current) / dtprev^order_current;
							% The p derivative of y at last time step
							fpnm1 = temp2 * nordsieck(:, order_current + 1) / dtprev^(order_current + 1);
							% The p-1 derivative of y at next time step
							fn = temp * nordsieck_next(:, order_current) / dtcur^order_current;
							% the p derivative of y at next time step
							fpn = temp2 * nordsieck_next(:, order_current + 1) / dtprev^(order_current + 1);

							% Estimated p + 1 derivative at last time step - can use this for order p error estimate
							backward_diff = (fn - fnm1 - dtcur * fpnm1) / dtcur^2;
							% Estimated p + 1 derivative at next time step - can use this for filling nordsieck vector
							forward_diff = (dtcur * fpn - fn + fnm1) / dtcur^2;

							% Estimated p + 2 derivative - can use this for order p + 1 error estimate
							final_diff = (forward_diff - backward_diff) / dtcur;

							% 24 from 4! from 4-point finite difference above
							err = errNorm.errEstimate(ycur, ynext, obj.C(order_current + 1) * final_diff * 24);
							[~, dtnext_candidate] = stepController.newStepSize(prevAccept, tspan, dtcur, err, order_current + 1);

							if dtnext_candidate > dtnext
								dtnext = dtnext_candidate;
								order_next = order_current + 1;
								nordsieck_next(:, order_current + 2) = forward_diff * 6 / temp2 / (order_current + 1) * dtcur ^ (order_current + 1);
							end
						end
					else
						% TODO setup with memory based time step controller
						prevAccept = false;
						% TODO: Allow factor to be choosen
						% TODO: Lighter penalty for if we are using an older Jacobian
						dtnext = 0.25 * dtcur;

						if order_current > 1
							order_next = order_current - 1;
						end
					end

					% Set new step to be in range
					dtnext = max(abs(dtmin), min(abs(dtmax), abs(dtnext))) * tdir;

					% Advance time
					if prevAccept
						% Kahan summation - adjust dt by adding the parts of previous dt
						% that we haven't been able to add yet
						dtadj = dtcur + dtbuf;
						tnext = tcur + dtadj;

						% Due to rounding error, tnext - tcur may not properly capture all of dtadj,
						% so record the parts that are missing to be added later
						dtbuf = dtadj - (tnext - tcur);
						break;
					end

					if abs(dtcur) == abs(dtmin)
						% TODO maybe we should return our progress up until now
						error("OneStepIntegrator:MinStep", "Step failed with h = MinStep")
					end

					stats.nFailed = stats.nFailed + 1;

					omega = dtnext / dtcur;
					nordsieck = nordsieck * sparse(1:(obj.MaxOrder+1), 1:(obj.MaxOrder+1), omega .^ (0:obj.MaxOrder));

					% If we increase order, ensure that the higher nordsieck elements are zeroed
					if order_next > order_current
						nordsieck(:, (order_current + 2) : (order_next + 1)) = zeros(length(y0), order_next - order_current);
					end

					dtcur = dtnext;
					order_current = order_next;
				end

				stats.nSteps = stats.nSteps + 1;

				% TODO: Add FullTrajectory and dense output

				if tcur * tdir >= tspan(tindex) * tdir
					tindex = tindex + 1;
				end

				if tspan(tindex) * tdir <= (tnext + dtnext) * tdir

					%integrate to/ End point
					%check if close enough with hmin
					if abs(tspan(tindex) - tnext) < 64 * eps(tnext)

						if multiTspan
							t(:, tindex) = tspan(tindex);
							y(:, tindex) = ynext;
						end

						if tcur < tspan(tindex)
							tindex = tindex + 1;
						end
					else

						% TODO - what if this is < dtmin
						dtnext = tspan(tindex) - tnext;
					end

				end
			end

			y(:, end) = ynext;
			t(end) = tspan(:, end);
		end

		function [err, stats] = timeStepErr(obj, y, ynew, ~, ErrNorm, order, delta, stats, opts)
			yerror = obj.C(order) / (obj.C(order) + 1) * obj.L(order, 1) * delta;

			if opts.StiffCorrectError
				[yerror, stats] = obj.NonLinearSolver.LinearSolver.solve(yerror, stats)
			end

			err = ErrNorm.errEstimate(y, ynew, yerror);
		end

		function [t, y, stats] = timeLoopFixed(obj, f, tspan, y0, opts)

			if opts.FullTrajectory
				y = zeros(length(y0), length(tspan));
				y(:, 1) = y0;
				t = tspan;
			else
				y = zeros(length(y0), 2);
				t = zeros(size(tspan,1),2);
			end

			t(:,1) = tspan(:,1);
			y(:, 1) = y0;

			%inital values
			ynext = y0;

			%Start stats
			stats = obj.intalizeStats;
			stats.nSteps = length(tspan);

			% nordsieck vector stored as a multivector, one column for each entry
			% TODO - get more values from startup procedure
			nordsieck = zeros(length(y0), obj.MaxOrder + 1);
			nordsieck(:, 1) = y0;
			nordsieck(:, 2) = (tspan(2) - tspan(1)) * f.F(tspan(1), y0);
			stats.nFevals = 1;

			% TODO - can get this from startup procedure
			current_order = 1;

			%Time Loop
			for i = 1:(length(tspan)-1)
				yi = ynext;
				tcur = tspan(i);
				dtnext = tspan(i+1) - tspan(i);

				if i > 1
					omega = dtnext / dtc;
				else
					omega = 1;
				end

				dtc = dtnext;

				Omega = sparse(1:(obj.MaxOrder+1), 1:(obj.MaxOrder+1), omega .^ (0:obj.MaxOrder));
				nordsieck = nordsieck * Omega;

				[ynext, nordsieck, ~, stats] = obj.timeStep(f, tcur, yi, dtc, nordsieck, true, current_order, stats);

				% Adaptive order strategy for fixed step - increase order every step until we hit max order
				% TODO - replace this with user-configured max order
				current_order = min(current_order + 1, obj.MaxOrder);

				if opts.FullTrajectory
					y(:, i + 1) = ynext;
					t(i + 1) = tspan(:, i + 1);
				end
			end

			y(:, end) = ynext;
			t(end) = tspan(:, end);
		end

		function [ynew, nordsieck, delta, stats, out_opts] = timeStep(obj, f, t, y, dt, nordsieck, prevAccept, order, stats)

			%% Predictor step - propagate nordsieck vector forward in time.
			% For explicit methods, the predicted y value (nordsieck(:, 1)) is exactly the next step value
			% Since nordsieck vector is a multivector, this is normal matrix multiplication from the left
			nordsieck(:, 1:(order + 1)) = nordsieck(:, 1:(order + 1)) * obj.Pascal(1:(order + 1), 1:(order + 1));

			l = obj.L(order, :);

			[~, stats] = obj.NonLinearSolver.preprocess(f, t, y, 1, dt * l(1), [], stats);

			%% Corrector step - correct explicit y value and higher derivatives
			% nordsieck(:,1) is a good choice of initial guess for the Nonlinear Solve since it is the next step of an explicit method
			% y0 must be nordsieck(:, 1), so any adjustment to the initial guess must be through x0, which is added to y0.
			x0 = zeros(size(y));
			sys_const = -l(1) * nordsieck(:, 2);
			[ydiff, solver_opts, stats] = obj.NonLinearSolver.solve(f, t + dt, dt, nordsieck(:, 1), x0, [], sys_const, 1, dt * l(1), [], stats);
			ynew = nordsieck(:, 1) + ydiff;

			nordsieck(:, 1) = ynew;

			fnew = f.F(t + dt, ynew);
			stats.nFevals = stats.nFevals + 1;

			delta = (dt * fnew - nordsieck(:, 2));
			nordsieck(:, 2:(order+1)) = nordsieck(:, 2:(order+1)) + delta * l(2:(order+1));

			if solver_opts.convergenceFailure == true
				out_opts.failure = true;
				return;
			end

			out_opts.failure = false;
		end

		function stats = intalizeStats(obj)

			stats.nFevals = 0;
			stats.nSteps = 0;
			stats.nFailed = 0;
			stats.nSmallSteps = 0;
			stats.nLinearSolves = 0;
			stats.nMassEvals = 0;
			stats.nJacobianEvals = 0;
			stats.nNonLinIterations = 0;
			stats.nPDTEval = 0;
			stats.nDecompositions = 0;
		end
	end
end