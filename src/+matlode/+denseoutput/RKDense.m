classdef RKDense < matlode.denseoutput.DenseOutput
    %RKDENSE Used for dense output computed from RK trees
    
    properties (SetAccess = immutable)
        bTheta
    end
    
    methods
        function obj = RKDense(b)
            obj.bTheta = @(theta) b * theta.^(1:size(b,2))';
        end
        
        function [denseY, fEvals] = denseOut(obj, ~, t, tneed, y0, ~, stages, dt)
            
            theta = (tneed - t) / dt;
            
            denseY = y0 + (stages * dt) * obj.bTheta(theta);
            fEvals = 0;
        end
    end
end

