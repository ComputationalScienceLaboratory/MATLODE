classdef Fixed < matlode.stepsizecontroller.StepSizeController
    %Fixed stepsize controller
    
    properties (Constant)
        Adaptive = false;
    end
    
    properties (SetAccess = immutable, GetAccess = private)
        Steps
    end
    
    methods
        
        function obj = Fixed(steps)
            obj = obj@matlode.stepsizecontroller.StepSizeController(1);
            obj.Steps = steps;
        end
        
        function [h0, f0, fevals] = startingStep(obj, ~, tspan, ~, ~, ~, ~, ~)
            h0 = (tspan(2) - tspan(1)) / obj.Steps;
            f0 = [];
            fevals = 0;
        end
        
        function [accept, h, tNew] = newStepSize(~, ~, ~, tspan, h, ~, ~, nSteps, ~, ~)
            accept = true;
            tNew = (nSteps) * h + tspan(1);
        end
        
    end
end

