classdef NordsieckDense < matlode.denseoutput.DenseOutput
	%NordsieckDense Used for dense output computed from nordsieck vector of multistep methods

	methods

		function [denseY, fEvals] = denseOut(obj, ~, t, tneed, y0, y1, dense_opts, dt)

			assert(tneed <= t + dt)

			fEvals = 0;
			theta = (tneed - t) / dt;
			dtf0 = dense_opts.nordsieck{1}(:, 2);
			dtf1 = dense_opts.nordsieck{2}(:, 2);
			denseY = (1 - theta) * y0 + theta * y1 + theta*(theta - 1) * ((1-2*theta) * (y1 - y0) + (theta - 1) * dtf0 + theta * dtf1);

			if dense_opts.order >= 4
				fp0 = dense_opts.nordsieck{1}(:, 3);
				fp1 = dense_opts.nordsieck{2}(:, 3);

				denseY = denseY + theta*theta*(theta - 1)*(theta-1)*(-3*(1 - 2*theta) * (y1 - y0) + (2 - 3*theta) * dtf0 + (1 - 3 * theta) * dtf1 + (1 - theta) * fp0 + theta * fp1);
			end
		end
	end
end

