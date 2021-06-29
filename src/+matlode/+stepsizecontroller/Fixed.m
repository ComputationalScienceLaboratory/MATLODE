classdef Fixed < matlode.stepsizecontroller.StepSizeController
    %Fixed stepsize controller
    
    properties (Constant)
        Adaptive = false;
    end
    
    properties (SetAccess = immutable, GetAccess = private)
        Steps
    end
    
    properties (SetAccess = immutable)
        History
    end
    
    methods
        
        function obj = Fixed(steps)
            obj.Steps = steps;
            obj.History = 1;
        end
        
        function [h0, f0] = startingStep(obj, f, tspan, y0, ~, ~, intialStep)
            if intialStep ~= 0
                warning('Ignoring initialStep option, Controller is Fixed');
            end
            
            f0 = f(tspan(1), y0);
            h0 = (tspan(2) - tspan(1)) / obj.Steps;
        end
        
        function [accept, h, tNew] = newStepSize(~, ~, ~, tspan, h, ~, ~, nSteps, ~, ~)
            accept = true;
            tNew = (nSteps) * h + tspan(1);
        end
        
    end
end

