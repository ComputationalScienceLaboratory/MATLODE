classdef ESDIRK_4_3_6 < matlode.rk.dirk.DIRK
	% ESDIRK 4(3)6
	% p = 4, s = 6, pe = 3
	
	methods
		function obj = ESDIRK_4_3_6(datatype)
            arguments
				datatype(1,1) string = 'double';
			end
            
            caster = @(x) matlode.util.CoefficentTransformers.transform(x,datatype);
            
            
            a = caster(join(['[[0, 0, 0, 0, 0, 0];',...
                              '[1/4, 1/4, 0, 0, 0, 0];',...
                              '[(1-sqrt(2))/8, (1-sqrt(2))/8, 1/4, 0, 0, 0];',...
                              '[(5 - 7*sqrt(2))/64, (5 - 7*sqrt(2))/64, (1 + sqrt(2))/32, 1/4, 0, 0];',...
                              '[(-13796 - 54539*sqrt(2))/125000, (-13796 - 54539*sqrt(2))/125000, (506605 + 132109*sqrt(2))/437500, 166*(-97 - 376*sqrt(2))/109375, 1/4, 0];',...
                              '[(1181 - 987*sqrt(2))/13782, (1181 - 987*sqrt(2))/13782, 47*(-267 - 1783*sqrt(2))/273343, -16*(-22922 + 3525*sqrt(2))/571953, -15625*(97 + 376*sqrt(2))/90749876, 1/4]]']));
             
            b =  caster('[(1181 - 987*sqrt(2))/13782, (1181 - 987*sqrt(2))/13782, 47*(-267 - 1783*sqrt(2))/273343, -16*(-22922 + 3525*sqrt(2))/571953, -15625*(97 + 376*sqrt(2))/90749876, 1/4]');
            
            bHat =  caster('[(263980 - 483197*sqrt(2))/4515000, (263980 - 483197*sqrt(2))/4515000, 483197*(1 + sqrt(2))/2257500, 293362/564375, -1/12, 10/43]');
            
            e = caster('[(367186009*sqrt(2))/10370955000 + 2352837/86424625, (367186009*sqrt(2))/10370955000 + 2352837/86424625, - (45894182153*sqrt(2))/88153117500 - 22915412153/88153117500, 13065461338/107598658125 - (18800*sqrt(2))/190651, 9070297/136124814 - (1468750*sqrt(2))/22687469, 3/172]');
            
            c = caster('[0, 1/2, (2 - sqrt(2))/4, 5/8, 26/25, 1]');
            
            order = 4;
            
            embbededOrder = 3;
            
            obj = obj@matlode.rk.dirk.DIRK(a, b, bHat, c, e, order, embbededOrder);
           
        end
        
        function sol = integrate(obj, f, tspan, y0, varargin)
            sol = integrate@matlode.rk.dirk.DIRK(obj, f, tspan, y0, 'StepSizeController', matlode.stepsizecontroller.StandardController, 'Dense', obj.DenseOut, varargin{:});
        end
	end
end

