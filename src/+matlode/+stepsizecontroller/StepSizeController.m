classdef (Abstract) StepSizeController < handle
    %Template for a step size controller
    
    properties (Abstract, Constant)
        Adaptive
    end
    
    properties (Abstract, SetAccess = immutable)
        History
    end
    
    methods (Abstract)
        %accept0, h0, and err0, respect the amount of history needed for a
        %controller
        [h0, f0] = startingStep(obj, f, tspan, y0, order, errFunc, intialStep);
        [accept, hNew, tNew] = newStepSize(obj, prevAccept, t, tspan, h, err, q, nSteps, nFailed);
    end
    
end

