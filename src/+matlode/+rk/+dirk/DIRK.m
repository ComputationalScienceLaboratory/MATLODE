classdef DIRK < matlode.rk.RungeKutta
	%DIRK Will support DIRK, ESDIRK, and SDIRK
	
    properties (Constant)
        PartitionMethod = false;
		PartitionNum = 1;
		MultirateMethod = false;
    end
    
    properties 
        NonLinearSolver
		StifflyAccurate
	end

	
	methods
        function obj = DIRK(a, b, bHat, c, e, order, embeddedOrder)
            obj = obj@matlode.rk.RungeKutta(a, b, bHat, c, e, order, embeddedOrder);
			
			obj.StifflyAccurate = all(obj.A(end,:) == obj.B(:)');
			if obj.FSAL
				
			end

        end
	end

	methods (Access = protected)
        function opts = matlodeSets(obj, p, varargin)
            
            %DIRK specific options
            p.addParameter('NonLinearSolver', matlode.nonlinearsolver.Chord());
            
            opts = matlodeSets@matlode.rk.RungeKutta(obj, p, varargin{:});

            if isempty(opts.NonLinearSolver)
                error('Please provide appropiate parameters and a non-linear solver')
            end
            obj.NonLinearSolver = opts.NonLinearSolver;
        end
        
        function [ynew, stages, stats, out_opts] = timeStep(obj, f, t, y, dt, stages, prevAccept, stats)
			if obj.FSAL && prevAccept
                stages(:, 1) = stages(:, end);
			end
			ynew = y;

			%Setup Nonlinear Solver
			[~, stats] = obj.NonLinearSolver.preprocess(f, t, y, [], stats);
            
			%TODO: Update PorMAss
			for i = obj.FsalStart:obj.StageNum
                g_const = 0;
                for j = 1:i-1
					if obj.A(i,j) ~= 0
						g_const = g_const + stages(:, j) .* (dt * obj.A(i, j));
					end
                end
                thc = t + dt .* obj.C(i);

				%Solve Nonlinear System
				[ydiff, solver_opts, stats] = obj.NonLinearSolver.solve(f, thc, y, zeros(length(y),1), g_const, 1, dt .* obj.A(i,i),  [], stats);
				if solver_opts.convergenceFailure == true
					out_opts.failure = true;
					return;
				end
				ynew = y + ydiff;

				stages(:,i) = f.F(thc, ynew);
			end

			if ~obj.StifflyAccurate
				ynew = y;
				for i = 1:obj.StageNum
					if obj.B(i) ~= 0
						ynew = ynew + (dt .* obj.B(i)) .* stages(:, i);
					end
				end
			end

			out_opts.failure = false;
            
		end
    end

end

