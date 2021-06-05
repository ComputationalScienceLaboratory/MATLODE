classdef ClassicRK4 < matlode.rk.erk.ERK
    methods
        function obj = ClassicRK4()
            
            %TODO
            %include reference
            
            a = [[0,   0,   0, 0];
                 [1/2, 0,   0, 0];
                 [0,   1/2, 0, 0];
                 [0,   0,   0, 1]];
             
            b = [1/6, 1/3, 1/3, 1/6];
            
            bHat = [];
            
            c = [0, 1/2, 1/2, 1];
            
            order = 4;
            
            embbededOrder = 0;
            
            obj = obj@matlode.rk.erk.ERK(a, b, bHat, c, order, embbededOrder);
        end
        
        function sol = integrate(obj, f, tspan, y0, varargin)
            sol = integrate@matlode.rk.erk.ERK(obj, f, tspan, y0, 'StepSizeController', matlode.stepsizecontroller.Fixed(1000), varargin{:});
        end
    end
end

