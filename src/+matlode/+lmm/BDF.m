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
		Pascal

		L % The matrix of methods coefficients. Each row is a different order method
		MaxOrder
	end

	methods
		function obj = BDF(datatype)
			arguments
				datatype(1,1) string = 'double';
			end

			obj = obj@matlode.Integrator(true, datatype);

			obj.MaxOrder = 5;
			obj.Pascal = zeros(obj.MaxOrder + 1, obj.MaxOrder + 1);
			obj.L = [...
				[1  , 1  , 0  , 0 , 0 , 0] / 1;...
				[2  , 3  , 1  , 0 , 0 , 0] / 3;...
				[6  , 11 , 6  , 1 , 0 , 0] / 11;...
				[24 , 50 , 35 , 10, 1 , 0] / 50;...
				[120, 274, 225, 85, 15, 1] / 274;...
				];

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

			opts = matlodeSets@matlode.Integrator(obj, p, varargin{:});

			if isempty(opts.NonLinearSolver)
				error('Please provide appropiate parameters and a non-linear solver')
			end
			obj.NonLinearSolver = opts.NonLinearSolver;
		end

		function [t, y, stats] = timeLoop(obj, f, tspan, y0, opts)
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
			% TODO - reuse this f call in the nonlinear solver
			nordsieck = zeros(length(y0), obj.MaxOrder + 1);
			nordsieck(:, 1) = y0;
			nordsieck(:, 2) = (tspan(2) - tspan(1)) * f.F(tspan(1), y0);
			stats.nFevals = 1;

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

				% Adaptive order strategy for fixed step - increase order every step until we hit max order
				currentorder = min(i, obj.MaxOrder);

				Omega = sparse(1:(obj.MaxOrder+1), 1:(obj.MaxOrder+1), omega .^ (0:obj.MaxOrder));
				nordsieck = nordsieck * Omega;

				[ynext, nordsieck, stats] = obj.timeStep(f, tcur, yi, dtc, nordsieck, true, currentorder, stats);

				if opts.FullTrajectory
					y(:, i + 1) = ynext;
					t(i + 1) = tspan(:, i + 1);
				end
			end

			y(:, end) = ynext;
			t(end) = tspan(:, end);
		end

		function [ynew, nordsieck, stats, out_opts] = timeStep(obj, f, t, y, dt, nordsieck, prevAccept, order, stats)

			% Predictor step - propagate nordsieck vector forward in time.
			% Since nordsieck vector is a multivector, this is normal matrix multiplication from the left
			nordsieck(:, 1:order + 1) = nordsieck(:, 1:order + 1) * obj.Pascal(1:order + 1, 1:order + 1);

			l = obj.L(order, :);

			[~, stats] = obj.NonLinearSolver.preprocess(f, t, y, 1, dt * l(1), [], stats);

			% TODO - using nordsieck(:, 1) as initial guess. Maybe try applying fixed point iteration first? Analyze the cost of doing so.
			% y0 must be nordsieck(:, 1), so any adjustment to the initial guess must be through x0, which is added to y0.
			% TODO - fn0 is passed in as [] in DIRK, only not if it's already precomputed in ESDIRK. Maybe can re-use this for fixed point iteration.
			x0 = zeros(size(y));
			sys_const = -l(1) * nordsieck(:, 2);
			[ydiff, solver_opts, stats] = obj.NonLinearSolver.solve(f, t + dt, dt, nordsieck(:, 1), x0, [], sys_const, 1, dt * l(1), [], stats);
			ynew = nordsieck(:, 1) + ydiff;

			nordsieck(:, 1) = ynew;

			fnew = f.F(t + dt, ynew);
			stats.nFevals = stats.nFevals + 1;

			nordsieck(:, 2:order+1) = nordsieck(:, 2:order+1) + (dt * fnew - nordsieck(:, 2)) * l(2:order+1);

			if solver_opts.convergenceFailure == true
				out_opts.failure = true;
				return;
			end
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