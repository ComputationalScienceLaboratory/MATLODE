classdef (Abstract) NormErr
    properties (SetAccess = immutable)
        RelTol
        AbsTol
    end
    
    methods (Abstract)
         errFunc = errEstimate(obj, y0, y1, yerror)
    end
    
    methods
        function obj = NormErr(absTol, relTol)
            obj.AbsTol = absTol;
            obj.RelTol = relTol;
        end
    end
end

