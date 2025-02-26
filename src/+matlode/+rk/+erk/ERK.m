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
        
        function [ynew, stages, stats, out_opts] = timeStep(obj, f, t, y, dt, stages, prevAccept, stats)
			
			if obj.FsalStart && prevAccept
                stages(:, 1) = stages(:, end);
			end
            
			%TODO: Update PorMass
            for i = obj.FsalStart:obj.StageNum
                ynew = y;
                for j = 1:i-1
					if obj.A(i,j) ~= 0
						ynew = ynew + stages(:, j) .* (dt .* obj.A(i, j));
					end
                end
                thc = t + dt .* obj.C(i);
                stages(:, i) = f.F(thc, ynew);
            end
            ynew = y;
			for i = 1:obj.StageNum
				if obj.B(i) ~= 0
					ynew = ynew + stages(:, i) .* (dt .* obj.B(i));
				end
			end
            
            stats.nFevals = stats.nFevals + double(obj.StageNum - obj.FsalStart + 1);

			out_opts.failure = false;
		end


        function [err, stats, out_opts] = timeStepErr(obj, ~, ~, y, ynew, dt, stages, ErrNorm, stats)
            
            
			%TODO: Update for Mass
            yerror = 0;
			for i = 1:obj.StageNum
				if abs(obj.E(i)) > eps
					yerror = yerror +  (dt .*  obj.E(i)) .* stages(:, i);
				end
			end
            err = ErrNorm.errEstimate(y, ynew, yerror);
			out_opts.failure = false;
        end
    end
end

