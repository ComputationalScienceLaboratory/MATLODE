classdef ERK < matlode.rk.RungeKutta
    %Explicit Runge Kutta class
    
    properties (Constant)
        PartitionMethod = false
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
        
        function [t, y, stats] = timeLoop(obj, f, tspan, y0, opts, varargin)
            
            numVars = length(y0);
            multiTspan = length(tspan) > 2;
            isDense = ~isempty(opts.Dense);
            overstep = false;
            startupErr = true;
            
            %Structure access optimizations
            errNormFunc = opts.ErrNorm;
            stepController = opts.StepSizeController;
            newStepFunc = @(prevAccept, tCur, tspan, hHist, err, q, nSteps, nFailed) opts.StepSizeController.newStepSize(prevAccept, tCur, tspan, hHist, err, q, nSteps, nFailed);
            
            
            k = zeros(numVars, obj.Stage);
            y = [];
            t = [];
            
            if multiTspan
                y = zeros(numVars, length(tspan));
                t = zeros(1, length(tspan));
            else
                y = zeros(numVars, opts.ChunkSize);
                t = zeros(1, opts.ChunkSize);
            end
            
            %inital values
            t(1,1) = tspan(1);
            tCur = tspan(1);
            tNext = tCur;
            ti = 2;
            
            tdir = sign(tspan(end) - tspan(1));
            hmax = opts.MaxStep * tdir;
            hmin = 64 * eps(tspan(1)) * tdir;
            
            y(:, 1) = y0;
            yNext = y0;
            
            %temporary stats
            nSteps = 1;
            nFailed = 0;
            nSmallSteps = 0;
            
            q = min(obj.Order, obj.EmbeddedOrder);
            
            %First step
            %allocates memory for first step
            [h0, f0, nFevals] = stepController.startingStep(f, tspan, y0, obj.Order, errNormFunc, hmin, hmax);
            
            hHist = ones(1, stepController.History) * h0;
            hN = hHist(1);
            err = zeros(1, length(hHist));
            
            accept = true;
            
            %%%%
            % Start of method specific intializing code
            %%%%
            
            fsalStart = uint32(obj.FSAL) + 1;
            fevalIterCounts = double(obj.Stage - fsalStart + 1);
            
            if obj.FSAL
                if nFevals == 0
                    k(:, end) = f(tspan(1), y0);
                    nFevals = nFevals + 1;
                else
                    k(:, end) = f0;
                end
            end
            
            %%%%
            % End of method specific intializing code
            %%%%
            
            %Time Loop
            while ti <= length(tspan)
                %Cycle history
                err = circshift(err, 1, 2);
                hHist = circshift(hHist, 1, 2);
                
                yCur = yNext;
                tCur = tNext;
                hC = hN;
                hHist(1) = hN;
                
                %Subject to change
                hmin = 64 * eps(tCur) * tdir;
                
                %%%%
                % Start of method specific before time loop code
                %%%%
                
                if obj.FSAL
                    k(:, 1) = k(:, end);
                end
                
                %%%%
                % End of method specific before time loop code
                %%%%
                
                %Accept Loop
                %Will keep looping until accepted step
                while true
                    
                    
                    %%%%
                    % Start of method specific time loop code
                    %%%%
                    
                    for i = fsalStart:obj.Stage
                        k(:, i) = f(tCur + hC * obj.C(i), yCur + k(:, 1:i-1) * (hC * obj.A(i, 1:i-1)'));
                    end

                    yNext = yCur + k * (hC * obj.B');

                    nFevals = nFevals + fevalIterCounts;
                    
                    prevAccept = accept;

                    %Perform error estimates
                    %Could maybe leave out of method specific zone?
                    if stepController.Adaptive
                        yEmbbeded =  k * (hC * obj.E');
                        err(:, 1) = errNormFunc([yNext, yCur], yEmbbeded);
                        
                        if startupErr
                            err(:, 2:end) = err(:, 1);
                        end
                        startupErr = false;
                    end
                    
                    %%%%
                    % End of method specific time loop code
                    %%%%
                    
                    %Find next Step
                    [accept, hN, tNext] = newStepFunc(prevAccept, tCur, tspan, hHist, err, q, nSteps, nFailed);
                    
                    %set new step to be in range
                    if stepController.Adaptive
                        hN = min(hmax * tdir, hN * tdir) * tdir;
                    end
                    
                    %Check if step is really small
                    if hN * tdir < hmin * tdir
                        if nSmallSteps == 0
                            warning('The step the integrator is taking extremely small, results may not be optimum')
                        end
                        %accept step since the step cannot get any smaller
                        nSmallSteps = nSmallSteps + 1;
                        hN = hmin;
                        tNext = tCur + hmin;
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
                
                nSteps = nSteps + 1;
                
                %check if controller failed to overstep for dense
                if isDense && overstep && tNext * tdir < tspan(ti) * tdir
                    overstep = false;
                end
                
                %Check for all memory allocation
                if ~multiTspan

                    %Allocate more memory if non-dense
                    if nSteps > length(t)
                        y = [y, zeros(numVars, opts.ChunkSize)];
                        t = [t, zeros(1, opts.ChunkSize)];
                    end

                    t(:, nSteps) = tNext;
                    y(:, nSteps) = yNext;
                end
                
                %check if condition holds
                %used for end checking,integrate to, and dense output
                if ~(tspan(ti) * tdir > (tNext + hN) * tdir)
                    %Dense Output
                    if isDense && multiTspan && ti ~= length(tspan)
                        
                        if overstep
                            while tNext * tdir > tspan(ti) * tdir && ti ~= length(tspan)
                                %Call Dense output function and record

                                t(:, ti) = tspan(ti);
                                [y(:, ti), tempFevals] = opts.Dense.denseOut(f, tCur, tspan(ti), yCur, yNext, k, hC);
                                nFevals = nFevals + tempFevals;

                                ti = ti + 1;
                            end
                            
                            if (tspan(ti) - tNext) * tdir < hmin * tdir
                                t(:, ti) = tspan(ti);
                                y(:, ti) = yNext;
                                ti = ti + 1;
                            end

                            overstep = false;
                        else
                            overstep = true;
                            
                            %check for tspan(end) to hit exactly
                            if stepController.Adaptive && ~(tspan(end) * tdir > (tNext + hN) * tdir)
                                hN = tspan(end) - tNext;
                            end
                        end
                    else
                        %integrate to/ End point
                        %check if close enough with hmin
                        if (tspan(ti) - tNext) * tdir < hmin * tdir
                            
                            if multiTspan
                                t(:, ti) = tspan(ti);
                                y(:, ti) = yNext;
                            end
                            ti = ti + 1;
                        else
                            if stepController.Adaptive
                                hN = tspan(ti) - tNext;
                            end
                            
                        end
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

