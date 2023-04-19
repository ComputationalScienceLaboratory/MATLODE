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
        
        function [ynew, stages, stats] = timeStep(obj, f, t, y, dt, stages, prevAccept, stats)
            persistent fsal s a b c
            if isempty(fsal)
                fsal = obj.FsalStart;
                s = obj.StageNum;
                a = obj.A;
                b = obj.B;
                c = obj.C;
            end
            
			if fsal && prevAccept
                stages(:, 1) = stages(:, end);
			end
            
			%TODO: Update PorMAss
			for i = fsal:s
                ynew = y;
                for j = 1:i-1
					if a(i,j) ~= 0
						ynew = ynew + stages(:, j) .* (dt * a(i, j));
					end
                end
                thc = t + dt .* c(i);

				%TODO: Setup NonLin Flags
				%Solve Nonlinear System
				[ynew, stages(:,i), ~, stats] = obj.NonLinearSolver.solve(f, thc, ynew, ynew, 1, dt .* a(i,i),  [], stats);
			end

			if ~obj.StifflyAccurate
				ynew = y;
				for i = 1:s
					if b(i) ~= 0
						ynew = ynew + stages(:, i) .* (dt .* b(i));
					end
				end
			end
            
		end
    end

end

