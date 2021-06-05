classdef DormandPrince < matlode.rk.erk.ERK
    methods
        function obj = DormandPrince()
            
            %TODO
            %include reference
            
            a = [[0,          0,           0,          0,        0,           0,     0];
                 [1/5,        0,           0,          0,        0,           0,     0];
                 [3/40,       9/40,        0,          0,        0,           0,     0];
                 [44/45,      -56/15,      32/9,       0,        0,           0,     0];
                 [19372/6561, -25360/2187, 64448/6561, -212/729, 0,           0,     0];
                 [9017/3168,  -355/33,     46732/5247, 49/176,   -5103/18656, 0,     0];
                 [35/384,     0,           500/1113,   125/192,  -2187/6784,  11/84, 0]];
             
            b = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0];
            
            bHat = [5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];
            
            c = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
            
            order = 5;
            
            embbededOrder = 4;
            
            obj = obj@matlode.rk.erk.ERK(a, b, bHat, c, order, embbededOrder);
        end
        
        function sol = integrate(obj, f, tspan, y0, varargin)
            sol = integrate@matlode.rk.erk.ERK(obj, f, tspan, y0, 'StepSizeController', matlode.stepsizecontroller.Fixed(1000), varargin{:});
        end
    end
end

