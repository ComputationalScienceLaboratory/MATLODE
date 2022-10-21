classdef ERK < matlode.rk.RungeKutta
    %Explicit Runge Kutta class
    
    properties (Constant)
        PartitionMethod = false;
		PartitionNum = 1;
		MultirateMethod = false;
    end

    methods
        function obj = ERK(a, b, bHat, c, e, order, embeddedOrder)
            obj = obj@matlode.rk.RungeKutta(a, b, bHat, c, e, order, embeddedOrder);
        end
    end
    
    methods (Access = protected)
        function opts = matlodeSets(obj, p, varargin)
            
            %ERK specific options
            
            opts = matlodeSets@matlode.rk.RungeKutta(obj, p, varargin{:});
        end
        
        function [ynew, stages, stats] = timeStep(obj, f, t, y, dt, stages, prevAccept, stats)
            persistent fsal s a b c fevalIterCounts
            if isempty(fsal)
                fsal = obj.FsalStart;
                s = obj.StageNum;
                a = obj.A;
                b = obj.B;
                c = obj.C;
                fevalIterCounts = double(obj.StageNum - obj.FsalStart + 1);
            end
            
			if fsal && prevAccept
                stages(:, 1) = stages(:, end);
			end
            
            for i = fsal:s
                ynew = y;
                for j = 1:i-1
					if a(i,j) ~= 0
						ynew = ynew + stages(:, j) .* (dt .* a(i, j));
					end
                end
                thc = t + dt .* c(i);
                stages(:, i) = f.F(thc, ynew);
            end
            ynew = y;
			for i = 1:s
				if b(i) ~= 0
					ynew = ynew + stages(:, i) .* (dt .* b(i));
				end
			end
            
            stats.nFevals = stats.nFevals + fevalIterCounts;
        end
        
        function [err, stats] = timeStepErr(obj, ~, ~, y, ynew, dt, stages, ErrNorm, stats)
            persistent e s
            if isempty(e)
                e = obj.E;
                s = obj.StageNum;
            end
            
            yerror = 0;
			for i = 1:s
				if e(i) ~= 0
					yerror = yerror + stages(:, i) .* (dt .* e(i));
				end
			end
            err = ErrNorm.errEstimate(y, ynew, yerror);
        end
    end
end

