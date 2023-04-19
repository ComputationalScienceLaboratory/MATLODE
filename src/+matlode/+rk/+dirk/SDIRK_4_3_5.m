classdef SDIRK_4_3_5 < matlode.rk.dirk.DIRK
	%SDIRK 4(3)5
	% p = 4, s = 5, pe = 3
	
	methods
		function obj = SDIRK_4_3_5(datatype)
            arguments
				datatype(1,1) string = 'double';
			end
            
            caster = @(x) matlode.util.CoefficentTransformers.transform(x,datatype);
            
            
            a = caster(join(['[[1/4, 0, 0, 0, 0];',...
                              '[13/20, 1/4, 0, 0, 0];',...
                              '[580/1287, -175/5148, 1/4, 0, 0];',...
                              '[12698/37375, -201/2990, 891/11500, 1/4, 0];',...
                              '[944/1365, -400/819, 99/35, -575/252, 1/4]]']));
             
            b =  caster('[944/1365, -400/819, 99/35, -575/252, 1/4]');
            
            bHat =  caster('[41911/60060, -83975/144144, 3393/1120, -27025/11088, 103/352]');
            
            e = caster('[-25/4004, 4525/48048, -45/224, 575/3696, -15/352]');
            
            c = caster('[1/4, 9/10, 2/3, 3/5, 1]');
            
            order = 4;
            
            embbededOrder = 3;
            
            obj = obj@matlode.rk.dirk.DIRK(a, b, bHat, c, e, order, embbededOrder);
           
        end
        
        function sol = integrate(obj, f, tspan, y0, varargin)
            sol = integrate@matlode.rk.dirk.DIRK(obj, f, tspan, y0, 'StepSizeController', matlode.stepsizecontroller.StandardController, 'Dense', obj.DenseOut, varargin{:});
        end
	end
end

