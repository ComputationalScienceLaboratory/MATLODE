classdef Direct < matlode.denseoutput.DenseOutput
    %This class is for custom denseoutputs that methods could require
    %The custum dense output cannot require f evaluations
    %Function must be in row format
    
    properties (SetAccess = immutable)
        kFuncs
    end
    
    methods
        function obj = Direct(kfuncs)
            obj.kFuncs = kfuncs;
        end
        
        function [denseY, fEvals] = denseOut(obj, ~, t, tneed, y0, ~, stages, dt)
            
            theta = (tneed - t) / h;
            
            denseY = y0 + (stages * dt) * obj.kFuncs(theta);
            fEvals = 0;
        end
    end
end

