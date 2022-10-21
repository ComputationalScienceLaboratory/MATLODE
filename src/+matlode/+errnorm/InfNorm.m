classdef InfNorm < matlode.errnorm.NormErr
   methods
       function obj = InfNorm(abstol, reltol)
           obj = obj@matlode.errnorm.NormErr(abstol, reltol);
       end
        
        function err = errEstimate(obj, y0, y1, yerror)
            sc = obj.AbsTol + max(abs(y0), abs(y1)) * obj.RelTol;
            err = norm(yerror ./ sc, inf);
        end
    end
end

