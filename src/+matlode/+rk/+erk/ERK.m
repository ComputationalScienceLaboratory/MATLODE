% classdef ERK < matlode.rk.RungeKutta
%     %Explicit Runge Kutta class
% 
%     properties (Constant)
%         PartitionMethod = false;
% 		PartitionNum = 1;
% 		MultirateMethod = false;
%     end
% 
%     properties
%         matrix_func
%     end
% 
%     methods
%         function obj = ERK(a, b, bHat, c, e, order, embeddedOrder)
%             obj = obj@matlode.rk.RungeKutta(a, b, bHat, c, e, order, embeddedOrder);
%         end
%     end
% 
%     methods (Access = protected)
%         function opts = matlodeSets(obj, p, varargin)
% 
%             %ERK specific options
%             p.addParameter('MatrixFunction', []);
% 
%             opts = matlodeSets@matlode.rk.RungeKutta(obj, p, varargin{:});
%             obj.matrix_func = opts.MatrixFunction;
%         end
% 
%         function [ynew, stages, stats, out_opts] = timeStep(obj, f, t, y, dt, stages, prevAccept, stats)
% 
% 			if obj.FsalStart && prevAccept
%                 stages(:, 1) = stages(:, end);
% 			end
% 
% 			%TODO: Update PorMass
%             for i = obj.FsalStart:obj.StageNum
%                 ynew = y;
%                 for j = 1:i-1
% 					if obj.A(i,j) ~= 0
% 						ynew = ynew + stages(:, j) .* (dt .* obj.A(i, j));
% 					end
%                 end
%                 thc = t + dt .* obj.C(i);
%                 stages(:, i) = f.F(thc, ynew);
%             end
%             ynew = y;
% 			for i = 1:obj.StageNum
% 				if obj.B(i) ~= 0
% 					ynew = ynew + stages(:, i) .* (dt .* obj.B(i));
% 				end
% 			end
% 
%             stats.nFevals = stats.nFevals + double(obj.StageNum - obj.FsalStart + 1);
% 
% 			out_opts.failure = false;
% 		end
% 
% 
%         function [err, stats, out_opts] = timeStepErr(obj, ~, ~, y, ynew, dt, stages, ErrNorm, stats)
% 
% 
% 			%TODO: Update for Mass
%             yerror = 0;
% 			for i = 1:obj.StageNum
% 				if abs(obj.E(i)) > eps
% 					yerror = yerror +  (dt .*  obj.E(i)) .* stages(:, i);
% 				end
% 			end
%             err = ErrNorm.errEstimate(y, ynew, yerror);
% 			out_opts.failure = false;
%         end
%     end
% end




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % CORRECTIONS % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
classdef ERK < matlode.rk.RungeKutta
    %Explicit Runge Kutta class

    properties (Constant)
        PartitionMethod = false;
		PartitionNum = 1;
		MultirateMethod = false;
    end

    properties
        matrix_func
        y_temp
    end

    methods
        function obj = ERK(a, b, bHat, c, e, order, embeddedOrder)
            obj = obj@matlode.rk.RungeKutta(a, b, bHat, c, e, order, embeddedOrder);
        end
    end

    methods (Access = protected)
        function opts = matlodeSets(obj, p, varargin)

            %ERK specific options
            p.addParameter('MatrixFunction', []);

            opts = matlodeSets@matlode.rk.RungeKutta(obj, p, varargin{:});
            obj.matrix_func = opts.MatrixFunction;
        end

        function [ynew, stages, stats, out_opts] = timeStep(obj, f, t, y, dt, stages, prevAccept, stats)

			if obj.FsalStart && prevAccept
                stages(:, 1) = stages(:, end);
			end

			%TODO: Update PorMass
            y_stage = zeros(length(y), obj.StageNum);
            for i = obj.FsalStart:obj.StageNum
                % ynew = y;
                y_stage(:, i) = y;
                for j = 1:i-1
					if obj.A(i,j) ~= 0
						% ynew = ynew + stages(:, j) .* (dt .* obj.A(i, j));
                        y_stage(:, i) = y_stage(:, i) + stages(:, j) .* (dt .* obj.A(i, j));
					end
                end
                thc = t + dt .* obj.C(i);
                stages(:, i) = f.F(thc, y_stage(:, i));
            end
            ynew = y;
            for i = 1:obj.StageNum
				if obj.B(i) ~= 0
					ynew = ynew + stages(:, i) .* (dt .* obj.B(i));
				end
            end

            % START OF CORRECTIONS
            tol = 1.e-12;

            %predicted final solution
            y_p = max(ynew, tol);
            % y_p = max(y_stage(:,length(obj.B)), tol);

            %first, truncate versions of the stages are computed
            y_trunc = y_stage;
            for i = 1:length(obj.B)
                y_trunc(:,i) = max(y_stage(:,i), tol);
            end

            %second, averaged system matrix is computed
            g_corr = zeros(length(y), length(y));
            for idx = 1:obj.StageNum
                g_corr = g_corr + obj.B(idx)*obj.matrix_func(y_trunc(:,idx), thc)*diag(ones(50,1)./y_p);
            end

            %corrected solution for y
            y_corr = (eye(length(g_corr)) - dt * g_corr)\y;
            % y_corr = max(y_corr, tol);

            obj.y_temp = ynew;
            ynew = y_corr;

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
            % err = ErrNorm.errEstimate(y, ynew, yerror);
            err = ErrNorm.errEstimate(y, obj.y_temp, yerror);
			out_opts.failure = false;
        end
    end
end


