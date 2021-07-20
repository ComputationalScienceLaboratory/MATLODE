classdef StandardNorm < matlode.errnorm.NormErr
    methods
        function obj = StandardNorm(abstol, reltol)
           obj = obj@matlode.errnorm.NormErr(abstol, reltol);
        end
        
        function err = errEstimate(obj, y, yE)
            persistent abstol reltol
            if isempty(abstol)
                abstol = obj.AbsTol;
                reltol = obj.RelTol;
            end
            sc = abstol + max(abs(y(:, 1)), abs(y(:, 2))) * reltol;
            err = sqrt(1/size(y,1)) * norm(yE ./ sc);
        end
    end
end

