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
            
            %Structure access optimizations
            errNormFunc = opts.ErrNorm;
            stepController = opts.StepSizeController;
            newStepFunc = @(prevAccept, tCur, tspan, hHist, err, q, nSteps, nFailed) opts.StepSizeController.newStepSize(prevAccept, tCur, tspan, hHist, err, q, nSteps, nFailed);
            
            
            
            k = zeros(numVars, obj.Stage);
            y = [];
            t = [];
            
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
            
            q = min(obj.Order, obj.EmbeddedOrder);
            
            %First step
            %allocates memory for first step
            [h0, f0, nFevals] = stepController.startingStep(f, tspan, y0, obj.Order, errNormFunc, hmin, hmax);
            
            hHist = ones(1, stepController.History) * h0;
            hN = hHist(1);
            err = ones(1, length(hHist));
            
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
                    
                    khc = k * hC;
                    yNext = yCur + khc * obj.B';

                    %Perform error estimates
                    %Could maybe leave out of method specific zone?
                    if stepController.Adaptive
                        yEmbbeded =  khc * obj.E';
                        err(:, 1) = errNormFunc([yNext, yCur], yEmbbeded);
                    end
                    
                    %%%%
                    % End of method specific time loop code
                    %%%%
                    
                    nFevals = nFevals + fevalIterCounts;
                    
                    prevAccept = accept;
                    
                    %Find next Step
                    [accept, hN, tNext] = newStepFunc(prevAccept, tCur, tspan, hHist, err, q, nSteps, nFailed);
                    
                    %Check if step is really small
                    if abs(hN) < abs(hmin)
                        if nSmallSteps == 0
                            warning('MATLODE:smallStepSize', 'The step the integrator is taking extremely small, results may not be optimal')
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
                
                %set new step to be in range
                if stepController.Adaptive
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

