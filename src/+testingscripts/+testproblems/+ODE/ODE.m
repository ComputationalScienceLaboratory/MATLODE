classdef (Abstract) ODE
    %ODE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = public)
        problem_name
        y_0
        t_0
    end

    methods (Abstract)
        dy = f(obj, t,y);
        df_y = jac(obj, t,y);
        df_t = f_t(obj, t,y);

        y = y_exact(obj, t);
    end

    
end

