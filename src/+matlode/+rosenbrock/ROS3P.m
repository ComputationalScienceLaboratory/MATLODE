classdef ROS3P < matlode.rosenbrock.Rosenbrock
    properties
       Dense 
    end
    
    methods
        function obj = ROS3P(datatype)
            if nargin == 0
                datatype = 'double';
            end
            
            caster = @(x) matlode.util.CoefficentTransformers.transform(x,datatype);
            
			gammadia = caster('[7.886751345948129e−01, 7.886751345948129e−01, 7.886751345948129e−01]');
            alphasum = caster('[0, 1, 0]');
            gammasum = caster('[7.886751345948129e−01, −2.113248654051871e−01, −1.077350269189626]');
            
            a = caster(join(['[[0, 0, 0];',...
                             '[1.267949192431123, 0, 0];',...
                             '[1.267949192431123, 0, 0]']));
            
            c = caster(join(['[[0, 0, 0];',...
                             '[−1.607695154586736, 0, 0];',...
                             '[−3.464101615137755, −1.732050807568877, 0]']));
            m = caster('[2, 5.773502691896258e−01, 4.226497308103742e−01]');
            me = caster('[2.113248654051871, 1,  4.226497308103742e−01]');
            
            order = 3;
            
            embbededOrder = 2;
            
            obj = obj@matlode.rosenbrock.Rosenbrock(gammadia, gammasum, alphasum, a, c, m, me, order, embbededOrder);
            
            obj.Dense = matlode.denseoutput.Linear(m);
        end
        
        function sol = integrate(obj, f, tspan, y0, varargin)
            sol = integrate@matlode.rosenbrock.Rosenbrock(obj, f, tspan, y0, 'StepSizeController', matlode.stepsizecontroller.StandardController, 'Dense', obj.DenseOut, varargin{:});
        end
    end
end

