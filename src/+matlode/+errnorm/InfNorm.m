classdef InfNorm < matlode.errnorm.NormErr
   methods
       function obj = InfNorm(abstol, reltol)
           obj = obj@matlode.errnorm.NormErr(abstol, reltol);
       end
        
        function err = errEstimate(obj, y, yE)
            persistent abstol reltol
            if isempty(abstol)
                abstol = obj.AbsTol;
                reltol = obj.RelTol;
            end
            sc = abstol + max(abs(y(:, 1)), abs(y(:, 2))) * reltol;
            err = norm(yE ./ sc, inf);
        end
    end
end

