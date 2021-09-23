classdef RK65SSP33 < matlode.rk.erk.ERK
    
    methods
        function obj = RK65SSP33(datatype)
            
            % A strong stability preserving RK method
            % Reference : Page 5.6, Table 4.12 a) https://core.ac.uk/download/pdf/56372752.pdf
            % Macdonald, Colin Barr. "Constructing high-order Runge-Kutta methods with embedded strong-stability-preserving pairs." 
            % PhD diss., Theses (Dept. of Mathematics)/Simon Fraser University, 2003.
            
            if nargin == 0
                datatype = 'double';
            end
            
            caster = @(x) matlode.util.CoefficentTransformers.transform(x,datatype);
            
            
            a = caster(join(['[[0,            0,            0,            0,            0,            0];',...
                             '[1,             0,            0,            0,            0,            0];',...
                             '[1/4,           1/4,          0,            0,            0,            0];',...
                             '[2046/15625,    -454/15625,   1533/15625,   0,            0,            0];',...
                             '[-739/5625,     511/5625,     -566/16875,   20/27,        0,            0];',...
                             '[11822/21875,   -6928/21875,  -4269/21875,  -4/7,         54/35,        0]]']));
             
            b =  caster('[1/6, 1/6, 2/3, 0, 0, 0]');
            
            bHat = caster('[1/24, 0, 0, 125/336, 27/56, 5/48]');
            
            e = caster('[1/8, 1/6, 2/3, -125/336, -27/56, -5/48]');
            
            c = caster('[0, 1, 1/2, 1/5, 2/3, 1]');
            
            order = 3;
            
            embbededOrder = 5;
            
            obj = obj@matlode.rk.erk.ERK(a, b, bHat, c, e, order, embbededOrder);
           
        end
        
        function sol = integrate(obj, f, tspan, y0, varargin)
            sol = integrate@matlode.rk.erk.ERK(obj, f, tspan, y0, 'StepSizeController', matlode.stepsizecontroller.StandardController, 'Dense', matlode.denseoutput.Linear(obj.B), varargin{:});
        end
    end
end

