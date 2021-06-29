classdef InfNorm < matlode.errnorm.NormErr
   methods (Static)
        function errFunc = errEstimate(AbsTol, RelTol)
            
            scFunc = @(y) (AbsTol + abs(y(:, 1)) * RelTol);
            errFunc = @(y, yE) norm(yE ./ scFunc(y), 'Inf');
            
        end
    end
end

