classdef (Abstract) NormErr
    methods (Abstract,Static)
         errFunc = errEstimate(AbsTol, RelTol)
    end
end

