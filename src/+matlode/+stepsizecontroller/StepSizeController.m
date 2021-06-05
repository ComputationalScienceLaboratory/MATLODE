classdef (Abstract) StepSizeController < handle
    %Template for a step size controller
    
    properties (Abstract, Constant)
        adaptive
    end
    
    methods (Abstract)
        
        [h0, f0] = startingStep(obj, f, tspan, y0, order, options);
        [accept, hNew, tNew] = newStepSize(obj, tCur, tspan, hCur, err, embeddedOrder, stepdlx, options);
    end
    
end

