classdef Rosenbrock < matlode.OneStepIntegrator
    %ROSENBROCK Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        PartitionMethod = false;
		PartitionNum = 1;
		MultirateMethod = false;
    end
    
    properties (SetAccess = private)
		GammaDia
        GammaSum
        A
        AlphaSum
        C
        M
        E
        Order
        StageNum
        EmbeddedOrder
        FSAL
    end

    properties (SetAccess = immutable, GetAccess = private)
        FsalStart
	end

	properties (SetAccess = protected)
        DenseOut
	end
    
    properties 
        LinearSolver
	end

	methods (Static)
		%% Coefficent Transformation Methods
		%Coefficent transformer
		function [gammadia, gammasum, alphasum, a, c, m, me, e] = RosCoefMethTrans(gamma, alpha, b, be)
			gammadia = diag(gamma);
			gammasum = sum(gamma, 2);
			alphasum = sum(alpha, 2);

			m = mrdivide(gamma, b);
			%TODO: Fix
			if isempty(be)
				me = [];
				e = [];
			else
				me = mrdivide(be, gamma);
				e = m - me;
			end

			eta = gamma * diag(1 ./ gammadia) + eye(length(gammadia));

			c = gamma;
			a = alpha;

			for i = 1:length(gamma)
				c(:,i) = gamma \ eta(:,i);
			end

			for i = 1:length(gamma)
				a(i,:) = mrdivide(alpha(i,:), gamma);
			end

		end

		%Coefficent transformer
		function [gamma, alpha, b, be] = RosCoefBaseTrans(gammadia, a, c, m, me)
			gamma = c;

			eta = eye(length(gammadia));
			xi = diag(gammadia) - c;
			for i = 1:length(c)
				gamma(:,i) = xi \ eta(:,i);
			end

			alpha = a * gamma;

			b = m * gamma;
			if isempty(me)
				be = [];
			else
				be = me * gamma;
			end
			

		end
	end
    
    methods

        function obj = Rosenbrock(gammadia, gammasum, alphasum, a, c, m, e, order, embeddedOrder)
            
            obj = obj@matlode.OneStepIntegrator(~isempty(e), class(gammadia));

			obj.GammaDia = gammadia;
			obj.GammaSum = gammasum;
			obj.AlphaSum = alphasum;

			obj.A = a;
			obj.C = c;
			obj.M = m;
			obj.E = e;
            
            obj.EmbeddedOrder = embeddedOrder;
            obj.Order = order;
            obj.StageNum = size(obj.M, 2);
            obj.FSAL = all(obj.A(end, :) == obj.M) && all(obj.A(1, :) == 0) && all(obj.C(1, :) == 0) && obj.AlphaSum(1) == 0;
            obj.FsalStart = uint32(obj.FSAL) + 1;
            obj.DenseOut = matlode.denseoutput.Linear(obj.M);
		end
        
	end
    
    methods (Access = protected)
        function opts = matlodeSets(obj, p, varargin)
            
            %Rosenbrock sepcific options
            p.addParameter('LinearSolver', matlode.linearsolver.MatrixLinearSolver());
            
            opts = matlodeSets@matlode.OneStepIntegrator(obj, p, varargin{:});
            
            if isempty(opts.LinearSolver)
                error('Please provide appropiate parameters and a linear solver')
            end
            obj.LinearSolver = opts.LinearSolver;
            
        end
        
		%% Time step
		function [ynew, stages, stats] = timeStep(obj, f, t, y, dt, stages, prevAccept, stats)
            
            if obj.FSAL && prevAccept
                stages(:, 1) = stages(:, end);
            end
            
			%Check if time derivative is avalible
			if ~isempty(f.PartialDerivativeTime)
            	dfdt_0 = f.PartialDerivativeTime(t, y);
				stats.nPDTEval = stats.nPDTEval + 1;
			end
            
            for i = obj.FsalStart:obj.StageNum
                ynew = y;
                for j = 1:(i-1)
                    if obj.A(i,j) ~= 0
                        ynew = ynew + obj.A(i,j) * stages(:,j);
                    end
                end

                ynew = f.F(t + obj.AlphaSum(i) * dt, ynew);

                for j = 1:(i-1)
                    if obj.C(i,j) ~= 0
                        ynew = ynew + (obj.C(i,j) / dt) * stages(:,j);
                    end
                end


				if ~isempty(f.PartialDerivativeTime) && obj.GammaSum(i) ~= 0
					ynew = ynew + obj.GammaSum(i) * dt * dfdt_0;
				end

				% Update Linear Solver
				if i == obj.FsalStart || obj.GammaDia(i) ~= obj.GammaDia(i-1)
					stats = obj.LinearSolver.preprocess(f, t, y, i==obj.FsalStart, 1/(dt * obj.GammaDia(i)), -1, stats);
				end

                [stages(:, i), stats] = obj.LinearSolver.solve(ynew, stats);
            end

            ynew = y;
            for i = 1:obj.StageNum
                if obj.M(i) ~= 0
                    ynew = ynew + obj.M(i) * stages(:,i);
                end
            end
            
            stats.nFevals = stats.nFevals + obj.StageNum - uint16(prevAccept);
        end
        
		%% Error Estimate
		function [err, stats] = timeStepErr(obj, ~, ~, y, ynew, ~, stages, ErrNorm, stats)

            y_error = 0;
            for i = 1:obj.StageNum
                if obj.E(i) ~= 0
                    y_error = y_error + obj.E(i) * stages(:,i);
                end
            end
            
            err = ErrNorm.errEstimate(y, ynew, y_error);
        end
        
        function [stages, stats] = timeLoopBeforeLoop(obj, f, f0, t0, y0, stats)
           stages = zeros(length(y0), obj.StageNum);
            
            if obj.FSAL
                if isempty(f0)
                    stages(:, end) = f.F(t0, y0);
                    stats.FEvals = stats.FEvals + 1;
                else
                    stages(:, end) = f0;
                end
            end
            
        end
        
        function [q] = timeLoopInit(obj)
            q = min(obj.Order, obj.EmbeddedOrder);
        end
    end
end

