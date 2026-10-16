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
			% obj.C = [1/2, 1/3, 1/4, 1/5, 1/6];

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

			%Start stats
			stats = obj.intalizeStats;
			stats.nSteps = 1;

			tlen = 0;
			if multiTspan
				y = zeros(numVars, length(tspan));
				t = zeros(1, length(tspan));
				stats.order = zeros(1, length(tspan));
			elseif opts.FullTrajectory
				y = zeros(numVars, opts.ChunkSize);
				t = zeros(1, opts.ChunkSize);
				stats.order = zeros(1, opts.ChunkSize);
				tlen = length(t);
			else
				y = zeros(length(y0), 2);
				t = tspan;
				stats.order = zeros(1, 2);
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
			delta = zeros(size(y0));


			%First step
			%allocates memory for first step
			% TODO - can get this from startup procedure
			order_next = 1;
			order_current = 1;
			[dt0, f0, stats.nFevals] = stepController.startingStep(f, tspan, y0, order_next, errNorm, dtmin, dtmax);
			dtnext = dt0;
			dtcur = dt0;

			% nordsieck vector stored as a multivector, one column for each entry
			% TODO - get more values from startup procedure
			nordsieck_next = zeros(length(y0), obj.MaxOrder + 1);
			nordsieck_next(:, 1) = y0;
			nordsieck_next(:, 2) = dtnext * f0;

			prevAccept = true;

			steps_since_order_change = 0;

			% time loop
			while tindex <= tspanlen
				if(stats.nSteps > opts.MaxNumSteps)
					warning("OneStepIntegrator:MaxNumSteps", "Integrator reached maximum steps before finishing integration at t = %f.", tcur)
					break
				end

				if order_next ~= order_current || dtnext ~= dtcur
					steps_since_order_change = 0;
				end

				ycur = ynext;
				tcur = tnext;
				dtprev = dtcur;
				dtcur = dtnext;
				omega = dtcur / dtprev;
				nordsieck = nordsieck_next * sparse(1:(obj.MaxOrder+1), 1:(obj.MaxOrder+1), omega .^ (0:obj.MaxOrder));
				order_current = order_next;
				delta_prev = delta;

				% Accept Loop
				% Will keep looping until accepted step
				while true
					[ynext, nordsieck_next, delta, stats, out_opts] = obj.timeStep(f, tcur, ycur, dtcur, nordsieck, true, order_current, stats);

					if out_opts.failure == false
						[err, stats] = timeStepErr(obj, ycur, ynext, dtcur, errNorm, order_current, delta, stats, opts);

						% Find next Step
						% TODO: Update to take stats
						[prevAccept, dt_opt] = stepController.newStepSize(prevAccept, tspan, dtcur, err, order_current);
						order_opt = order_current;

						% Check the step size if we decrease order
						if order_current > 1
							% Radhakrishnan & Hindmarsh "Description and use of LSODE" p. 34
							% The nordsieck vector already stores an estimate of y_p, which is the residual for order p - 1
							yerr = obj.C(order_current - 1) * factorial(order_current) * nordsieck_next(:, order_current + 1);

							% TODO: Consider adding LSODE's different safety factors for different orders?
							%       Maybe we can do this by having different step controllers for each order.
							err = errNorm.errEstimate(ycur, ynext, yerr);
							[~, dtnext_candidate] = stepController.newStepSize(prevAccept, tspan, dtcur, err, order_current - 1);

							if dtnext_candidate > dt_opt
								dt_opt = dtnext_candidate;
								order_opt = order_current - 1;
							end
						end
					else
						% TODO setup with memory based time step controller
						prevAccept = false;
						% TODO: Allow factor to be choosen
						% TODO: Lighter penalty for if we are using an older Jacobian
						dtnext = 0.25 * dtcur;

						% TODO: Record when we decrease order because of nonlinear fail v.s. error test
						if order_current > 1
							steps_since_order_change = 0;
							order_next = order_current - 1;
						end
					end

					% Check to adapt order for the next step
					if prevAccept && steps_since_order_change > order_current
						if order_current < obj.MaxOrder
							% Radhakrishnan & Hindmarsh "Description and use of LSODE" p. 34
							% Use finite difference on delta from previous step to get higher derivative
							yerr = obj.C(order_current + 1) * factorial(order_current) * obj.L(order_current, order_current + 1) * (delta - delta_prev);

							err = errNorm.errEstimate(ycur, ynext, yerr);
							[~, dtnext_candidate] = stepController.newStepSize(prevAccept, tspan, dtcur, err, order_current + 1);

							% fprintf("Considering increase in order: %e v.s. %e\n", dtnext_candidate, dtnext);
							if dtnext_candidate > dt_opt
								dt_opt = dtnext_candidate;
								order_opt = order_current + 1;

								% Radhakrishnan & Hindmarsh "Description and use of LSODE" p. 36
								nordsieck_next(:, order_current + 2) = obj.L(order_current, order_current + 1) * delta / (order_current + 1);
							end
						end
					end

					if dt_opt < dtcur || dt_opt > 1.5 * dtcur
						order_next = order_opt;
						dtnext = dt_opt;
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

					if order_next ~= order_current || dtnext ~= dtcur
						steps_since_order_change = 0;
					end

					dtcur = dtnext;
					order_current = order_next;
				end

				stats.nSteps = stats.nSteps + 1;
				steps_since_order_change = steps_since_order_change + 1;

				if ~multiTspan && opts.FullTrajectory

					%Allocate more memory if non-dense
					if stats.nSteps > tlen
						y = [y, zeros(numVars, opts.ChunkSize)];
						t = [t, zeros(1, opts.ChunkSize)];
						stats.order = [stats.order, zeros(1, opts.ChunkSize)];
						tlen = length(t);
					end

					t(:, stats.nSteps) = tnext;
					y(:, stats.nSteps) = ynext;
					stats.order(:, stats.nSteps - 1) = order_current;

					while tindex <= tspanlen && tnext * tdir > tspan(tindex) * tdir
						tindex = tindex + 1;
					end
				else
					%Dense output
					while tindex <= tspanlen && tnext * tdir > tspan(tindex) * tdir

						%incase lucky and land on point, can be cheaper
						%than dense output if dense requires fevals
						if abs(tspan(tindex) - tnext) < 64 * eps(tcur)
							t(:, tindex) = tspan(tindex);
							y(:, tindex) = ynext;
						else
							%Call Dense output function and record
							t(:, tindex) = tspan(tindex);

							dense_opts.order = order_current;
							dense_opts.nordsieck = {nordsieck, nordsieck_next};
							dense = matlode.denseoutput.NordsieckDense;
							[y(:, tindex), tempFevals] = dense.denseOut(f, tcur, tspan(tindex), ycur, ynext, dense_opts, dtcur);
							stats.nFevals = stats.nFevals + tempFevals;
						end
						tindex = tindex + 1;
					end
				end
			end

			if ~multiTspan && opts.FullTrajectory
				t = t(:, 1:stats.nSteps);
				y = y(:, 1:stats.nSteps);
				stats.order = stats.order(:, 1:stats.nSteps);
				stats.order(end) = order_next;
			end
		end

		function [err, stats] = timeStepErr(obj, y, ynew, ~, ErrNorm, order, delta, stats, opts)
			% yerror = obj.C(order) / (obj.C(order) + 1) * obj.L(order, 1) * delta;
			yerror = obj.C(order) * factorial(order) * obj.L(order, order + 1) * delta;

			if opts.StiffCorrectError
				[yerror, stats] = obj.NonLinearSolver.LinearSolver.solve(yerror, stats);
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
			stats.order = [];
		end
	end
end