classdef DIRK < matlode.rk.RungeKutta
	%DIRK Will support DIRK, ESDIRK, and SDIRK
	
    properties (Constant)
        PartitionMethod = false;
		PartitionNum = 1;
		MultirateMethod = false;

		%Using Z_Transformation varies. At times false can be faster even
		%with the more Feval. The nonlinear system seems better
		%conditioned. However, z transformation works best for DAEs. SDIRK
		%can work for M(t) DAEs only under Z-transformation
		Z_Transformation = true;
    end
    
    properties 
        NonLinearSolver
	end

    properties (SetAccess = protected)
		Gamma
        D
		Theta
		ET
		Alpha
		StifflyAccurate
		ESDIRK
		SinglyImplicit

		falseYDiff
	end

	
	methods
        function obj = DIRK(a, b, bHat, c, e, order, embeddedOrder)
            obj = obj@matlode.rk.RungeKutta(a, b, bHat, c, e, order, embeddedOrder);

			adiag = diag(obj.A);
			
			obj.StifflyAccurate = all(obj.A(end,:) == obj.B(:)');
			%TODO: Change to coefficents provides to prevent accuracy loss.
			%Similar to the Rosenbrock methods

			%Note is only for SDIRK and ESDIRK methods. FIRK and other
			%will require more care. Firk requries full inverses
			%ESDIRK only utilizes the sub matrices inverse
			obj.ESDIRK = all(obj.A(1, :) == 0);
			if (obj.ESDIRK)
				obj.Gamma = obj.A(2,2);

				obj.SinglyImplicit = all(abs(adiag(2:end) - obj.Gamma) < eps);

				%Note the first stage is treated differently in the z
				%transformation. Adding the effects to the intial stage are
				%nesscary
				q = -(obj.A(2:end, 2:end) \ obj.A(2:end, 1));

				obj.D = zeros(1,obj.StageNum);
				obj.D(1) = obj.B(1) + obj.B(1,2:end)*q;
				if obj.StifflyAccurate
					obj.D(1,end) = 1;
				else
					obj.D(1,2:end) = (obj.A(2:end, 2:end)' \ obj.B(2:end)')';
				end
	
				obj.ET = zeros(1, obj.StageNum);
				if obj.Adaptive
					obj.ET(1) = obj.E(1) + obj.E(1,2:end)*q;
					obj.ET(1,2:end) = (obj.A(2:end, 2:end)' \ (obj.E(2:end))')';
				else
					obj.ET = [];
				end
	
				obj.Theta = zeros(obj.StageNum, obj.StageNum);
				obj.Theta(2:end, 2:end) = obj.A(2:end, 2:end) \ tril(obj.A(2:end, 2:end), -1);
				obj.Theta(2:end, 1) = obj.A(2:end,1) + tril(obj.A(2:end, 2:end), -1)*q;

			else
				obj.Gamma = obj.A(1,1);

				obj.SinglyImplicit = all(abs(adiag - obj.Gamma) < eps );

				if obj.StifflyAccurate
					obj.D = zeros(1,obj.StageNum);
					obj.D(1,end) = 1;
				else
					obj.D = (obj.A' \ obj.B')';
				end
	
				if obj.Adaptive
					obj.ET = (obj.A' \ (obj.E)')';
				else
					obj.ET = obj.E;
				end

				obj.Theta = obj.A \ tril(obj.A, -1);
			end

			if ~obj.SinglyImplicit
				error('General Dirk is not supported. Only SDIRK and ESDIRK.');
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
			if isa(f.Mass,'function_handle') && (~obj.Z_Transformation || obj.ESDIRK)
				error('Not Supported method for time dependent Mass Matrix M(t)')
			end

			if obj.FSAL && prevAccept
				if obj.Z_Transformation
					stages(:, 1) = f.F(t,y);
					
					stats.nFevals = stats.nFevals + 1;
				else
					stages(:, 1) = stages(:, end);
				end
			end
			ynew = y;

			%Setup Nonlinear Solver
			%TODO:Currently assuming diagonals are the same. fix
			[~, stats] = obj.NonLinearSolver.preprocess(f, t, y, 1, dt .* obj.Gamma, [], stats);
            
			if obj.Z_Transformation
				ydiff = zeros(length(y),1);
				fd0 = [];
				for i = obj.FsalStart:obj.StageNum

					g_const = zeros(length(y),1);
                	for j = 2:i-1
						if abs(obj.Theta(i,j)) > eps
							g_const = g_const + (obj.Theta(i, j)) .* stages(:, j) ;
						end
                	end
                	thc = t + obj.C(i) .* dt ;

					if ~obj.ESDIRK && i > 1 && abs(obj.Theta(i,1)) > eps
						g_const = g_const + obj.Theta(i,1) .* stages(:,1);
					end

					[stats] = obj.NonLinearSolver.LinearSolver.computeMass(f, thc, y, stats);
					g_const = (obj.NonLinearSolver.LinearSolver.mass * g_const);

					if obj.ESDIRK && i > 1 && abs(obj.Theta(i,1)) > eps
						g_const = g_const +  dt * obj.Theta(i,1) .* stages(:,1);
					end
	
					%Solve Nonlinear System
					[stages(:,i), solver_opts, stats] = obj.NonLinearSolver.solve(f, thc, dt, y, ydiff, fd0, g_const, 1, dt .* obj.Gamma,  [], stats);
					ydiff = stages(:,i);
					obj.falseYDiff(:,i) = ydiff;
					if solver_opts.convergenceFailure == true
						out_opts.failure = true;
						return;
					end
				end
				
				if obj.ESDIRK
					if abs(obj.D(1)) > eps
						ynew = ynew + dt * (obj.D(i)) .* stages(:, 1);
					end
				end
	
				if ~obj.StifflyAccurate
					esdirkstart = uint32(obj.ESDIRK) + 1;
					for i = esdirkstart:obj.StageNum
						if abs(obj.D(i)) > eps
							ynew = ynew + (obj.D(i)) .* stages(:, i);
						end
					end
				else
					ynew = ynew + stages(:,end);
				end
			else
				obj.falseYDiff = zeros(length(y),obj.StageNum);
				ydiff = zeros(length(y),1);
				if obj.ESDIRK
					fd0 = stages(:,1);
				else
					fd0 = [];
				end

				for i = obj.FsalStart:obj.StageNum
                	g_const = zeros(length(y),1);
                	for j = 1:i-1
						if obj.A(i,j) ~= 0
							g_const = g_const + (dt * obj.A(i, j)) .* stages(:, j) ;
						end
                	end
                	thc = t + obj.C(i) .* dt;
	
					%Solve Nonlinear System
					[ydiff, solver_opts, stats] = obj.NonLinearSolver.solve(f, thc, dt, y, ydiff, fd0, g_const, 1, dt .* obj.Gamma,  [], stats);
					if solver_opts.convergenceFailure == true
						out_opts.failure = true;
						return;
					end
					ynew = y + ydiff;
					obj.falseYDiff(:,i) = ydiff;
	
					stages(:,i) = f.F(thc, ynew);
					fd0 = stages(:,i);
				end
				stats.nFevals = stats.nFevals + obj.StageNum;
	
				if ~obj.StifflyAccurate
					ynew = y;
					for i = 1:obj.StageNum
						if obj.B(i) ~= 0
							ynew = ynew + (dt .* obj.B(i)) .* stages(:, i);
						end
					end
				end

			end

			out_opts.failure = false;
            
		end


        function [err, stats, out_opts] = timeStepErr(obj, f, t, y, ynew, dt, stages, ErrNorm, stats)
            yerror = 0;
			if obj.Z_Transformation
				if obj.ESDIRK
					if abs(obj.ET(1)) > eps
						yerror = dt*(obj.ET(1)) .* stages(:, 1);
					end
				end
				esdirkstart = uint32(obj.ESDIRK) + 1;
				for i = esdirkstart:obj.StageNum
					if abs(obj.ET(i)) > eps
						yerror = yerror +  (obj.ET(i)) .* stages(:, i);
					end
				end
				
			else
				% for i = 1:obj.StageNum
				% 	if abs(obj.E(i)) > eps
				% 		yerror = yerror +  (dt .*  obj.E(i)) .* stages(:, i);
				% 	end
				% end

				if obj.ESDIRK
					if abs(obj.ET(1)) > eps
						yerror = dt*(obj.ET(1)) .* stages(:, 1);
					end
				end
				esdirkstart = uint32(obj.ESDIRK) + 1;
				
				for i = esdirkstart:obj.StageNum
					if abs(obj.ET(i)) > eps
						yerror = yerror + (obj.ET(i)) .* obj.falseYDiff(:, i);
					end
				end
				% [stats] = obj.NonLinearSolver.LinearSolver.computeMass(f, t + dt, ynew, stats);
				% yerror = (obj.NonLinearSolver.LinearSolver.mass * yerror);
				% [yerror, stats] = obj.NonLinearSolver.LinearSolver.solve(yerror, stats);
			end

            err = ErrNorm.errEstimate(y, ynew, yerror);
			out_opts.failure = false;
        end
    end

end

