classdef (Abstract) OneStepIntegrator < matlode.Integrator
    %One Step Integrator template
    
    methods (Access = protected)
        function obj = OneStepIntegrator(varargin)
            obj = obj@matlode.Integrator(varargin{:});
        end
        
        function opts = matlodeSets(obj, p, varargin)
            
            %OneStepIntegrator Specific options
            
            opts = matlodeSets@matlode.Integrator(obj, p, varargin{:});
            
        end
        
		%% Time Loop
        function [t, y, stats] = timeLoop(obj, f, tspan, y0, opts)
            
            numVars = length(y0);
            multiTspan = length(tspan) > 2;
            isDense = ~isempty(opts.Dense);
            
            
            %Structure access optimizations
            errNorm = opts.ErrNorm;
            stepController = opts.StepSizeController;
            
            tlen = 0;
            
			% TODO: Switch to variable and not default
			if multiTspan
                y = zeros(numVars, length(tspan));
                t = zeros(1, length(tspan));
			elseif opts.FullTrajectory
                y = zeros(numVars, opts.ChunkSize);
                t = zeros(1, opts.ChunkSize);
                tlen = length(t);
			else
                y = zeros(length(y0), 2);
                t = tspan;
			end
			t(:,1) = tspan(:,1);
			y(:, 1) = y0;
            
            %inital values
            tcur = tspan(1);
            tnext = tcur;
            dtbuf = 0; % Kahan summation buffer to remember small parts of dt which haven't yet been added to t
            tindex = 2;
            
            tspanlen = length(tspan);
            
            tdir = sign(tspan(end) - tspan(1));
            dtmax = min([abs(opts.MaxStep), abs(tspan(end) - tspan(1))]) * tdir;
            dtmin = opts.MinStep * tdir;
            
            ynext = y0;
            
			%Create Stats Object
			stats = obj.intalizeStats();
			stats.nSteps = 1;
            
            [q] = obj.timeLoopInit;
            
            %First step
            %allocates memory for first step
            [dt0, f0, stats.nFevals] = stepController.startingStep(f, tspan, y0, q, errNorm, dtmin, dtmax);
            
            dthist = ones(1, stepController.History) * dt0;
            dtnext = dthist(1);
            err = ones(1, length(dthist));
            
            prevAccept = true;
            
            [stages, stats] = obj.timeLoopBeforeLoop(f, f0, tspan(1), y0, stats);
            
            %Time Loop
            while tindex <= tspanlen
                if(stats.nSteps > opts.MaxNumSteps)
                    warning("OneStepIntegrator:MaxNumSteps", "Integrator reached maximum steps before finishing integration at t = %f.", tcur)
                    break
                end

                %Cycle history
                err = circshift(err, 1, 2);
                dthist = circshift(dthist, 1, 2);
                
                ycur = ynext;
                tcur = tnext;
                dtcur = dtnext;
                dthist(1) = dtnext;
                
                %Accept Loop
                %Will keep looping until accepted step
                while true
                    
                    [ynext, stages, stats, out_opts] = timeStep(obj, f, tcur, ycur, dtcur, stages, prevAccept, stats);
					%Check if integration failed with the given time step
					if out_opts.failure == false
                    	[err(:, 1), stats, out_opts] = timeStepErr(obj, f, tcur, ycur, ynext, dtcur, stages, errNorm, stats);
                    	
                    	%Find next Step
						%TODO: Update to take stats
                    	[prevAccept, dtnext] = stepController.newStepSize(prevAccept, tspan, dthist, err, q);
					else
						%TODO setup with memory based time step controller
						prevAccept = false;
						%TODO: Allow factor to be choosen
						dtnext = 0.1 * dtcur;
					end

                	% Verify we don't go below MinStep.
                	if abs(dtnext) < abs(dtmin)
                    	dtnext = dtmin;
                	end
                	
                	%check step acception
                    if prevAccept
                        % Kahan summation - adjust dt by adding the parts of previous dt
                        % that we haven't been able to add yet
                        dtadj = dtcur + dtbuf;
                        tnext = tcur + dtadj;

                        % Due to rounding error, tnext - tcur may not properly capture all of dtadj,
                        % so record the parts that are missing to be added later
                        dtbuf = dtadj - (tnext - tcur);
                        break;
                    else
                        if dtcur == dtmin
                            % TODO maybe we should return our progress up until now
                            error("OneStepIntegrator:MinStep", "Step failed with h = MinStep")
                        end
                    end
                	
                	stats.nFailed = stats.nFailed + 1;
                	dtcur = dtnext;
                	dthist(1) = dtnext;

                end
                
                %set new step to be in range
                dtnext = min(abs(dtmax), abs(dtnext)) * tdir;
                
                stats.nSteps = stats.nSteps + 1;
                
                %Check for all memory allocation
                if ~multiTspan && opts.FullTrajectory

                    %Allocate more memory if non-dense
                    if stats.nSteps > tlen
                        y = [y, zeros(numVars, opts.ChunkSize)];
                        t = [t, zeros(1, opts.ChunkSize)];
                        tlen = length(t);
                    end

                    t(:, stats.nSteps) = tnext;
                    y(:, stats.nSteps) = ynext;
                elseif isDense
                    %Dense output
                    while tindex ~= tspanlen && tnext * tdir > tspan(tindex) * tdir 
                        
                        %incase lucky and land on point, can be cheaper
                        %than dense output if dense requires fevals
                        if abs(tspan(tindex) - tnext) < 64 * eps(tcur)
                            t(:, tindex) = tspan(tindex);
                            y(:, tindex) = ynext;
                        else
                            %Call Dense output function and record
                            t(:, tindex) = tspan(tindex);
                            [y(:, tindex), tempFevals] = opts.Dense.denseOut(f, tcur, tspan(tindex), ycur, ynext, stages, dtcur);
                            stats.nFevals = stats.nFevals + tempFevals;
                        end
                        tindex = tindex + 1;
                    end
                end
                
                %check if condition holds
                %used for end checking and integrate to
                if tspan(tindex) * tdir <= (tnext + dtnext) * tdir
                    
                    %integrate to/ End point
                    %check if close enough with hmin
                    if abs(tspan(tindex) - tnext) < 64 * eps(tnext)

                        if multiTspan
                            t(:, tindex) = tspan(tindex);
                            y(:, tindex) = ynext;
                        end
                        tindex = tindex + 1;
					else
                        
                        % TODO - what if this is < dtmin
                        dtnext = tspan(tindex) - tnext;
                    end
                    
                end
                
                
            end
            
            %End of Timeloop work
            %Truncate extra memory
            if ~multiTspan && opts.FullTrajectory
                t = t(:, 1:stats.nSteps);
                y = y(:, 1:stats.nSteps);
            else
                t(:, end) = tspan(end);
                y(:, end) = ynext;
            end
            
            %Create Stats
            stats.nSteps = stats.nSteps - 1;
            
		end

		%% Fixed Time Loop
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


            [stages, stats] = obj.timeLoopBeforeLoop(f, [], tspan(1), y0, stats);
            
            %Time Loop
			for i = 1:(length(tspan)-1)
                yi = ynext;
                tcur = tspan(i);
                dtc = tspan(i+1) - tspan(i);
                
                [ynext, stages, stats] = timeStep(obj, f, tcur, yi, dtc, stages, true, stats);

				if opts.FullTrajectory
					y(:, i + 1) = ynext;
                    t(i + 1) = tspan(:, i + 1);
				end
			end

			y(:, end) = ynext;
            t(end) = tspan(:, end);
		end

		%Statistics intialization
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

