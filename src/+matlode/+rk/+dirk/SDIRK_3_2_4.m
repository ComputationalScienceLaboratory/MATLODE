classdef SDIRK_3_2_4 < matlode.rk.dirk.DIRK
	%SDIRK 3(2)4
	% p = 3, s = 4, pe = 2
	
	methods
		function obj = SDIRK_3_2_4(datatype)
            arguments
				datatype(1,1) string = 'double';
			end
            
            caster = @(x) matlode.util.CoefficentTransformers.transform(x,datatype);
            
            
            a = caster(join(['[[9/40, 0, 0, 0];'...
                              '[163/520, 9/40, 0, 0];',...
                              '[-6481433/8838675, 87795409/70709400, 9/40, 0];',...
                              '[4032/9943, 6929/15485, -723/9272, 9/40] ]']));
             
            b =  caster('[4032/9943, 6929/15485, -723/9272, 9/40]');
            
            bHat =  caster('[20/51, 64477140871/138472716300, -1303583701/18463028840, 1034014989/4858691800]');
            
            e = caster('[6772/507093, -410011317313/22571052756900, -2075562131/281561189810, 29595333/2429345900]');
            
            c = caster('[9/40, 7/13, 11/15, 1]');
            
            order = 3;
            
            embbededOrder = 2;
            
            obj = obj@matlode.rk.dirk.DIRK(a, b, bHat, c, e, order, embbededOrder);
           
        end
        
        function sol = integrate(obj, f, tspan, y0, varargin)
            sol = integrate@matlode.rk.dirk.DIRK(obj, f, tspan, y0, 'StepSizeController', matlode.stepsizecontroller.StandardController, 'Dense', obj.DenseOut, varargin{:});
        end
	end
end

