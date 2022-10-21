classdef StandardNorm < matlode.errnorm.NormErr
    methods
        function obj = StandardNorm(abstol, reltol)
           obj = obj@matlode.errnorm.NormErr(abstol, reltol);
        end
        
        function err = errEstimate(obj, y0, y1, yerror)
            sc = obj.AbsTol + max(abs(y0), abs(y1)) * obj.RelTol;
            err = sqrt(1/size(y0,1)) * norm(yerror ./ sc);
        end
    end
end

