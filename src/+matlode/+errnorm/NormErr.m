classdef (Abstract) NormErr
    properties (SetAccess = immutable)
        RelTol
        AbsTol
    end
    
    methods (Abstract)
         errFunc = errEstimate(y, yE)
    end
    
    methods
        function obj = NormErr(absTol, relTol)
            obj.AbsTol = absTol;
            obj.RelTol = relTol;
        end
    end
end

