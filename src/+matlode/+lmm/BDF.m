classdef BDF < matlode.Integrator
	properties (Constant)
		PartitionMethod = false;
		PartitionNum = 1;
		MultirateMethod = false;
	end

	properties (SetAccess = protected)
		Pascal

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

			opts = matlodeSets@matlode.Integrator(obj, p, varargin{:});

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


			[nordsieck, stats] = obj.timeLoopBeforeLoop(f, [], tspan(1), y0, stats);

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

				if abs(omega - 1) > 1e-14
					error("BDF currently only works for fixed step integration")
				end

				dtc = dtnext;

				[ynext, nordsieck, stats] = obj.timeStep(f, tcur, yi, dtc, nordsieck, true, stats);

				% TODO - Omega matrix. Because of check above we know omega is always approximately 1, so Omega must be identity.
				% So we can skip multiplying here. To get rid of the check above and allow variable time steps, we must add Omega matrix multiplication
				% of Nordsieck vector.

				if opts.FullTrajectory
					y(:, i + 1) = ynext;
					t(i + 1) = tspan(:, i + 1);
				end
			end

			y(:, end) = ynext;
			t(end) = tspan(:, end);
		end

		function [ynew, nordsieck, stats, out_opts] = timeStep(obj, f, t, y, dt, nordsieck, prevAccept, stats)

			% Predictor step - propagate nordsieck vector forward in time.
			% Since nordsieck vector is a multivector, this is normal matrix multiplication from the left
			nordsieck = nordsieck * obj.Pascal;


		end

		function [nordsieck, stats] = timeLoopBeforeLoop(obj, f, f0, t0, y0, stats)

			% nordsieck vector stored as a multivector, one column for each entry
			nordsieck = zeros(length(y0), obj.MaxOrder + 1);

		end

		function [q] = timeLoopInit(obj)

			q = obj.MaxOrder;
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