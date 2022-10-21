classdef BogackiShampine < matlode.rk.erk.ERK
%Method: Dormand PRince
% p = 3 s = 4 pe = 2
%Reference: 
    
    methods

        function obj = BogackiShampine(datatype)
            arguments
				datatype(1,1) string = 'double';
			end
            
            caster = @(x) matlode.util.CoefficentTransformers.transform(x,datatype);
            
            
            a = caster(join(['[[0, 0, 0, 0];',...
                              '[1/2, 0, 0, 0];',...
                              '[0, 3/4,	0, 0];',...
                              '[2/9, 1/3, 4/9, 0]]']));
             
            b =  caster('[2/9, 1/3, 4/9, 0]');
            
            bHat = caster('[7/24, 1/4, 1/3, 1/8]');
            
            e = caster('[-5/72, 1/12, 1/9, -1/8]');
            
            c = caster('[0, 1/2, 3/4, 1]');
            
            order = 3;
            
            embbededOrder = 2;
            
            obj = obj@matlode.rk.erk.ERK(a, b, bHat, c, e, order, embbededOrder);
           
        end
        
        function sol = integrate(obj, f, tspan, y0, varargin)
            sol = integrate@matlode.rk.erk.ERK(obj, f, tspan, y0, 'StepSizeController', matlode.stepsizecontroller.StandardController, 'Dense', obj.DenseOut, varargin{:});
        end
    end
end

