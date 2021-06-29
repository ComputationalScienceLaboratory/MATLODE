classdef H321 < matlode.stepsizecontroller.SoberlandController
    %TODO Reference NASA Paper and Soberland paper
    %This controller works and provides a smooth solution
    
    methods
        function obj = H321()
            a = [-1/3, -1/18, 5/18];
            b = [5/6, 1/6];
            
            AQfunc = @(q, A) (A / q);
            
            %This is stupid, it ignores the first parameter for some reason
            obj = obj@matlode.stepsizecontroller.SoberlandController(0, 'A', a, 'B', b, 'AQfunc', AQfunc);
        end
    end
end

