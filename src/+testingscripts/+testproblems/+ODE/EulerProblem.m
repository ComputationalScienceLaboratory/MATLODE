classdef EulerProblem < testingscripts.testproblems.ODE.ODE
    % The one dimensional Euler problem with lambda parameter
    
    properties
        lambda
    end
    
    methods
        function obj = EulerProblem(lambda)
            obj.problem_name = 'Euler Problem';
            obj.lambda = lambda;
            obj.y_0 = 1;
            obj.t_0 = 0;
        end

        function dy = f(obj, ~,y)
            dy = obj.lambda * y;
        end

        function df_y = jac(obj, ~,~)
            df_y = obj.lambda;
        end

        function df_t = f_t(~, ~,~)
            df_t = 0;
        end

        function y = y_exact(obj, t)
            y = exp(obj.lambda * t);
        end
    end
end

