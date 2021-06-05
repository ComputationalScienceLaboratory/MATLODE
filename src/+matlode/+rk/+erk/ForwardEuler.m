classdef ForwardEuler < matlode.rk.erk.ERK
    methods
        function obj = ForwardEuler()
            %TODO
            %include reference
            
            a = [[0]];
             
            b = [1];
            
            bHat = [];
            
            c = [0];
            
            order = 1;
            
            embbededOrder = 0;
            
            obj = obj@matlode.rk.erk.ERK(a, b, bHat, c, order, embbededOrder);
        end
        
        function sol = integrate(obj, f, tspan, y0, varargin)
            sol = integrate@matlode.rk.erk.ERK(obj, f, tspan, y0, 'StepSizeController', matlode.stepsizecontroller.Fixed(1000), varargin{:});
        end
    end
end

