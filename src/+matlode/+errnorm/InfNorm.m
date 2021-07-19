classdef InfNorm < matlode.errnorm.NormErr
   methods (Static)
        function errFunc = errEstimate(AbsTol, RelTol)
            
            scFunc = @(y) (AbsTol + max(abs(y(:, 1)), abs(y(:, 2))) * RelTol);
            errFunc = @(y, yE) norm(yE ./ scFunc(y), 'Inf');
            
        end
    end
end

