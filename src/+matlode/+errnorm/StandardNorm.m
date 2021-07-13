classdef StandardNorm < matlode.errnorm.NormErr
    methods (Static)
        function errFunc = errEstimate(AbsTol, RelTol)
            
            scFunc = @(y) (AbsTol + max(abs(y(:, 1)), abs(y(:, 2))) * RelTol);
            errFunc = @(y, yE) sqrt(1/size(y,1)) * norm(yE ./ scFunc(y));
            
        end
    end
end

