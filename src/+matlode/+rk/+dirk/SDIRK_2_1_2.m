classdef SDIRK_2_1_2 < matlode.rk.dirk.DIRK
	%SDIRK 2(1)2
	% p = 2, s = 2, pe = 1
	
	methods
		function obj = SDIRK_2_1_2(datatype)
            arguments
				datatype(1,1) string = 'double';
			end
            
            caster = @(x) matlode.util.CoefficentTransformers.transform(x,datatype);
            
            
            a = caster(join(['[[1 - (1/sqrt(2)), 0];',...
                              '[1/sqrt(2), 1 - (1/sqrt(2))]]']));
             
            b =  caster('[1/sqrt(2), 1 - (1/sqrt(2))]');
            
            bHat =  caster('[3/5, 2/5]');
            
            e = caster('[1/sqrt(2) - 3/5, 3/5 - 1/sqrt(2)]');
            
            c = caster('[1 - (1/sqrt(2)), 1]');
            
            order = 2;
            
            embbededOrder = 1;
            
            obj = obj@matlode.rk.dirk.DIRK(a, b, bHat, c, e, order, embbededOrder);
           
        end
        
        function sol = integrate(obj, f, tspan, y0, varargin)
            sol = integrate@matlode.rk.dirk.DIRK(obj, f, tspan, y0, 'StepSizeController', matlode.stepsizecontroller.StandardController, 'Dense', obj.DenseOut, varargin{:});
        end
	end
end

