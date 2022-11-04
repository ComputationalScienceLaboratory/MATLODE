classdef OrderChecking < handle
	%ORDERCHECKING Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
		Integrator
	end
	
	methods
		function obj = OrderChecking(int)
			obj.Integrator = int;
		end
		
		function [error, statsCumlative] = getErrorWRTExactValue(obj, model, t0, tf, y0, ytrue, stepamount, opts)
			
            error = zeros(1,length(stepamount));
			statsCumlative = cell(1, length(stepamount));
			

			for i = 1:length(error)
				steps = t0:((tf - t0) / stepamount(i)):tf;
				sol = obj.Integrator.integrateFixed(model, steps, y0, opts);
				statsCumlative{i} = sol.stats;
				error(i) = norm(sol.y(end) - ytrue) / norm(ytrue);
			end
		end

		function [f, polyError] = plot_error_single(obj, stepamount, error, fig_num)
            f = figure(fig_num);

            %% For Total Error
            loglog(stepamount, error)
            
            title(['Total Error of: ', class(obj.Integrator)])
            xlabel('Time step')
            ylabel('Relative Error')

            polyError = nan_fix(error);
            
            legend(['p = ', num2str(polyError(1))])

            function [p_cn] = nan_fix(error)
                %step amount becomes block variable
    
                %shift arrays to remove nan for poly fit for coefficents
                nanshift_t = rmmissing(error);
                err_size = length(error);
                
                %gives the linear polynomials of the errors
                p_cn = polyfit(log(stepamount((err_size + 1) - length(nanshift_t):err_size)), log(nanshift_t), 1);
            end
        end
	end
end

