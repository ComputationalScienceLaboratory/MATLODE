classdef (Abstract) StepSizeController < handle
    %Template for a step size controller
    
    properties (Abstract, Constant)
        Adaptive
    end
    
    properties (SetAccess = immutable)
        History
    end
    
    methods (Access = protected)
        function obj = StepSizeController(hist, varargin)
            
            p = inputParser;
            p.KeepUnmatched = false;
            
            %Gives better error checking
            p.parse(varargin{:});
            
           obj.History = hist; 
        end
    end
    
    methods (Abstract)
        %accept0, h0, and err0, respect the amount of history needed for a
        %controller
        [h0, f0, fevals] = startingStep(obj, f, tspan, y0, order, errFunc, minStep, maxStep);
        [accept, hNew, tNew] = newStepSize(obj, prevAccept, t, tspan, h, err, q, nSteps, nFailed);
    end
    
end

