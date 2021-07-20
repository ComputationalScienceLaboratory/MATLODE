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
        
        function [t, y, stats] = timeLoop(obj, f, tspan, y0, opts)
            
            numVars = length(y0);
            multiTspan = length(tspan) > 2;
            isDense = ~isempty(opts.Dense);
            
            
            %Structure access optimizations
            errNorm = opts.ErrNorm;
            stepController = opts.StepSizeController;
            adap = stepController.Adaptive;
            
            tlen = 0;
            
            if multiTspan
                y = zeros(numVars, length(tspan));
                t = zeros(1, length(tspan));
            else
                y = zeros(numVars, opts.ChunkSize);
                t = zeros(1, opts.ChunkSize);
                tlen = length(t);
            end
            
            %inital values
            t(1,1) = tspan(1);
            tCur = tspan(1);
            tNext = tCur;
            ti = 2;
            tl = 2;
            if isDense
               ti = length(tspan); 
            else
                tl = length(tspan);
            end
            
            tspanlen = length(tspan);
            
            tdir = sign(tspan(end) - tspan(1));
            hmax = min([abs(opts.MaxStep), abs(tspan(end) - tspan(1))]) * tdir;
            hmin = 64 * eps(tspan(1)) * tdir;
            
            y(:, 1) = y0;
            yNext = y0;
            
            %temporary stats
            nSteps = 1;
            nFailed = 0;
            nSmallSteps = 0;
            
            [q] = obj.timeLoopInit;
            
            %First step
            %allocates memory for first step
            [h0, f0, nFevals] = stepController.startingStep(f, tspan, y0, q, errNorm, hmin, hmax);
            
            hHist = ones(1, stepController.History) * h0;
            hN = hHist(1);
            err = ones(1, length(hHist));
            
            accept = true;
            
            [k, fevals, fevalIterCounts] = obj.timeLoopBeforeLoop(f0, tspan(1), y0);
            
            nFevals = nFevals + fevals;
            
            %Time Loop
            while ti <= tspanlen
                %Cycle history
                err = circshift(err, 1, 2);
                hHist = circshift(hHist, 1, 2);
                
                yCur = yNext;
                tCur = tNext;
                hC = hN;
                hHist(1) = hN;
                
                %Subject to change
                hmin = 64 * eps(tCur) * tdir;
                
                %Accept Loop
                %Will keep looping until accepted step
                while true
                    
                    prevAccept = accept;
                    
                    if adap
                        [k, yNext, err(:, 1)] = timeStepErr(obj, f, tCur, yCur, hC, k, errNorm, prevAccept);
                    else
                        [k, yNext, err(:, 1)] = timeStep(obj, f, tCur, yCur, hC, k, errNorm, prevAccept);
                    end
                    
                    nFevals = nFevals + fevalIterCounts;
                    
                    %Find next Step
                    [accept, hN, tNext] = stepController.newStepSize(prevAccept, tCur, tspan, hHist, err, q, nSteps, nFailed);
                    
                    %Check if step is really small
                    if abs(hN) < abs(hmin)
                        if nSmallSteps == 0
                            warning('The step the integrator is taking extremely small, results may not be optimal')
                        end
                        %accept step since the step cannot get any smaller
                        nSmallSteps = nSmallSteps + 1;
                        hN = hmin;
                        tNext = tCur + hmin;
                        accept = true;
                        break;
                    end
                    
                    %check step acception
                    if accept
                        break;
                    end
                    
                    nFailed = nFailed + 1;
                    hC = hN;
                    hHist(1) = hN;
                    
                end
                
                %set new step to be in range
                if adap
                    hN = min(abs(hmax), abs(hN)) * tdir;
                end
                
                nSteps = nSteps + 1;
                
                %Check for all memory allocation
                if ~multiTspan

                    %Allocate more memory if non-dense
                    if nSteps > tlen
                        y = [y, zeros(numVars, opts.ChunkSize)];
                        t = [t, zeros(1, opts.ChunkSize)];
                        tlen = length(t);
                    end

                    t(:, nSteps) = tNext;
                    y(:, nSteps) = yNext;
                elseif isDense
                    %Dense output
                    while tl ~= tspanlen && tNext * tdir > tspan(tl) * tdir 
                        
                        %incase lucky and land on point, can be cheaper
                        %than dense output if dense requires fevals
                        if abs(tspan(tl) - tNext) < abs(hmin)
                            t(:, tl) = tspan(tl);
                            y(:, tl) = yNext;
                        else
                            %Call Dense output function and record
                            t(:, tl) = tspan(tl);
                            [y(:, tl), tempFevals] = opts.Dense.denseOut(f, tCur, tspan(tl), yCur, yNext, k, hC);
                            nFevals = nFevals + tempFevals;
                        end
                        tl = tl + 1;
                    end
                end
                
                %check if condition holds
                %used for end checking and integrate to
                if tspan(ti) * tdir <= (tNext + hN) * tdir
                    
                    %integrate to/ End point
                    %check if close enough with hmin
                    if abs(tspan(ti) - tNext) < abs(hmin)

                        if multiTspan
                            t(:, ti) = tspan(ti);
                            y(:, ti) = yNext;
                        end
                        ti = ti + 1;
                    elseif stepController.Adaptive
                        
                        hN = tspan(ti) - tNext;
                    end
                    
                end
                
                
            end
            
            %End of Timeloop work
            %Truncate extra memory
            if ~multiTspan
                t = t(:, 1:nSteps);
                y = y(:, 1:nSteps);
            else
                t(:, end) = tspan(end);
                y(:, end) = yNext;
            end
            
            %Create Stats
            stats.nSteps = nSteps - 1;
            stats.nFailed = nFailed;
            stats.nFevals = nFevals;
            stats.nSmallSteps = nSmallSteps;
            
        end
        
    end
end

