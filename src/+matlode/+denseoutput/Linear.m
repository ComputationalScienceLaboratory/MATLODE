classdef Linear < matlode.denseoutput.DenseOutput
    
    properties (SetAccess = immutable)
       b 
    end
    
    methods
        function obj = Linear(b)
            obj.b = b;
        end
        
        function [denseY, fEvals] = denseOut(obj, ~, t, tNeeded, yC, ~, k, h)
            
            theta = (tNeeded - t) / h;
            denseY = yC + (k * h) * (obj.b' * theta) ;
            fEvals = 0;
        end
    end
end

