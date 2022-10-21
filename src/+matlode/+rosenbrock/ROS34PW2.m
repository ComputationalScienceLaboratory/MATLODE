classdef ROS34PW2 < matlode.rosenbrock.Rosenbrock
%Method: ROS34PW2
% Rosenbrock-w method
% p = 3 s = 4 pe = 2
%Reference: 
    
    methods
        function obj = ROS34PW2(datatype)
            arguments
				datatype(1,1) string = 'double';
			end
            
            caster = @(x) matlode.util.CoefficentTransformers.transform(x,datatype);
			
            gammadia = caster('[0.435866521508459, 0.435866521508459, 0.435866521508459, 0.435866521508459]');
            alphasum = caster('[0, 0.871733043016918, 0.731579957788852, 1]');
            gammasum = caster('[0.435866521508459, -0.435866521508459, -0.413333376233886, -5.55111512312578e-17]');
            
            a = caster(join(['[[0, 0, 0, 0];',...
                             '[2, 0, 0, 0];',...
                             '[1.41921731745576, -0.25923221167297, 0, 0];',...
                             '[4.18476048231916, -0.285192017355496, 2.29428036027904, 0]]']));
            
            c = caster(join(['[[0, 0, 0, 0];',...
                             '[-4.58856072055809, 0, 0, 0];',...
                             '[-4.18476048231916, 0.285192017355496, 0, 0];',...
                             '[-6.36817920012836, -6.79562094446684, 2.87009860433106, 0]]']));
            m = caster('[4.18476048231916, -0.285192017355496, 2.29428036027904, 1]');
            e = caster('[0.277749947647967, -1.403239895176, 1.77263012766755, 0.5]');
            
            order = 3;
            
            embbededOrder = 2;
            
            obj = obj@matlode.rosenbrock.Rosenbrock(gammadia, gammasum, alphasum, a, c, m, e, order, embbededOrder);
            
        end
        
        function sol = integrate(obj, f, tspan, y0, varargin)
            sol = integrate@matlode.rosenbrock.Rosenbrock(obj, f, tspan, y0, 'StepSizeController', matlode.stepsizecontroller.StandardController, 'Dense', obj.DenseOut, varargin{:});
        end
    end
end

