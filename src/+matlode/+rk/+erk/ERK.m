classdef ERK < matlode.rk.RungeKutta
    %Explicit Runge Kutta class
    
    properties (SetAccess = immutable)
        a
        b
        bHat
        c
        stage
        order
        embeddedOrder
        adaptable
        FSAL
    end
    
    methods
        
        
        
        function obj = ERK(a, b, bHat, c, order, embeddedOrder)
            
            obj.a = a;
            obj.b = b;
            obj.c = c;
            obj.bHat = bHat;
            obj.embeddedOrder = embeddedOrder;
            obj.order = order;
            obj.stage = length(b);
            obj.FSAL = all(a(end, :) == b) && all(a(1, :) == 0);
            obj.adaptable = ~isempty(bHat);
        end
        
        function sol = integrate(obj, f, tspan, y0, varargin)
            
            opts = obj.matlodeSets(varargin{:});
            
            numVars = length(y0);
            t = tspan(1);
            k = zeros(numVars, obj.stage);
            sol.y = y0;
            
            sol.stats.nsteps = 0;
            sol.stats.nfailed = 0;
            sol.stats.nfevals = 0;
            
            fsalStart = uint32(obj.FSAL) + 1;
            fevalIterCounts = obj.stage - fsalStart + 1;
            
            
            [h, k(:,1)] = opts.StepSizeController.startingStep(f, tspan, y0, obj.order, opts);
            
            
            while t < tspan(end)
                
                for i = fsalStart:obj.stage
                    k(:, i) = f(t + h * obj.c(i), sol.y + h * (k(:, 1:i-1) * obj.a(i, 1:i-1)'));
                end
                
                sol.stats.nfevals = sol.stats.nfevals + fevalIterCounts;
                
                yCurrent = sol.y + h * (k * obj.b');
                err = 0;
                
                
                if opts.StepSizeController.adaptive
                    yEmbbeded = sol.y + h * (k * obj.bHat);
                    err = opts.ErrNorm(yCurrent, yEmbbeded, opts);
                end
                
                [accept, h, t] = opts.StepSizeController.newStepSize(t, tspan, h, err, obj.embeddedOrder, sol.stats, opts);
                
                if accept
                    sol.y = yCurrent;
                    sol.stats.nsteps = sol.stats.nsteps + 1;
                    
                    
                    if obj.FSAL
                        k(:, 1) = k(:, end);
                    end
                else
                    
                   sol.stats.nfailed = sol.stats.nfailed + 1;
                end
            end
            
            
        end
        
        function opts = matlodeSets(obj, varargin)
            p = inputParser;
            
            %ERK specific options
            
            opts = obj.matlodeSetsHelper(p, varargin{:});
        end
    end
    
end

