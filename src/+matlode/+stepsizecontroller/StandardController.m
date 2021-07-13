classdef StandardController < matlode.stepsizecontroller.StartingController
    %Standard Error Controller as described in Solving ODES I book
    %Seperate from Soberland to help with performance
    
    
    methods
        function obj = StandardController(varargin)
            
            obj = obj@matlode.stepsizecontroller.StartingController(1, varargin{:});
            
        end
        
        function [accept, hNew, tNew] = newStepSize(obj, ~, t, ~, h, err, q, ~, ~)
            accept = err <= 1;
            
            if accept
                tNew = t + h;
            else
                tNew = t;
            end
            
            hNew = h * min(obj.FacMax, max(obj.FacMin, obj.Fac * err^(-1 / (q + 1))));
            
        end
    end
end

