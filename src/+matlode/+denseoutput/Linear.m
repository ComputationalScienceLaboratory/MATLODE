classdef Linear < matlode.denseoutput.DenseOutput
    
    properties (SetAccess = immutable)
       b 
    end
    
    methods
        function obj = Linear(b)
            obj.b = b;
        end
        
        function [denseY, fEvals] = denseOut(obj, ~, t, tneed, y0, ~, stages, dt)
            
            theta = (tneed - t) / dt;
            denseY = y0 + (stages * dt) * (obj.b' * theta) ;
            fEvals = 0;
        end
    end
end

