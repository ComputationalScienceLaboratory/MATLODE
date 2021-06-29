classdef StandardController < matlode.stepsizecontroller.StartingController
    %Standard Error Controller as described in Solving ODES I book
    %Seperate from Soberland to help with performance(TODO test if true)
    
    properties (Constant)
        Adaptive = true;
    end
    
    properties (SetAccess = immutable)
        Fac
        FacMin
        FacMax
        History
    end
    
    methods
        function obj = StandardController(varargin)
            
            obj.History = 1;
            p = inputParser;
            
            p.addParameter('Fac', 0.95);
            p.addParameter('FacMin', 0.1);
            p.addParameter('FacMax', 2);
            
            p.parse(varargin{:});
            
            opts = p.Results;
            
            obj.Fac = opts.Fac;
            obj.FacMin = opts.FacMin;
            obj.FacMax = opts.FacMax;
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

