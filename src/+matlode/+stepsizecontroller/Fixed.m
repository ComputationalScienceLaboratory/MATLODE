classdef Fixed < matlode.stepsizecontroller.StepSizeController
    %Fixed stepsize controller
    
    properties (Constant)
        adaptive = false;
    end
    
    properties (SetAccess = immutable, GetAccess = private)
        steps
    end
    
    methods
        
        function obj = Fixed(steps)
            obj.steps = steps;
        end
        
        function [h0, f0] = startingStep(obj, f, tspan, y0, ~, options)
            if options.InitialStep ~= 0
                warning('Ignoring initialStep option, Controller is Fixed');
            end
            
            f0 = f(tspan(1), y0);
            h0 = (tspan(2) - tspan(1)) / obj.steps;
        end
        
        function [accept, hCur, tNew] = newStepSize(~, ~, tspan, hCur, ~, ~, stepdlx, ~)
            accept = true;
            tNew = (stepdlx.nsteps + 1) * hCur + tspan(1);
        end
        
    end
end

