classdef LinBackEuler < matlode.rosenbrock.Rosenbrock
%Method: Linearly Implicit Euler
% p = 1 s = 1 pe = 0
%Reference: 
    
    methods
        function obj = LinBackEuler(datatype)
            arguments
				datatype(1,1) string = 'double';
			end
            
            caster = @(x) matlode.util.CoefficentTransformers.transform(x,datatype);
            
            alpha = caster('[0]');
            gamma = caster('[1]');
            b = caster('[1]');
            be = [];

			[gammadia, gammasum, alphasum, a, c, m, me, e] = RosCoefMethTrans(gamma, alpha, b, be);
            
            order = 1;
            
            embbededOrder = 0;
            
            obj = obj@matlode.rosenbrock.Rosenbrock(gammadia, gammasum, alphasum, a, c, m, e, order, embbededOrder);
            
        end
        
        function sol = integrate(obj, f, tspan, y0, varargin)
            sol = integrate@matlode.rosenbrock.Rosenbrock(obj, f, tspan, y0, 'Dense', obj.DenseOut, varargin{:});
        end
    end
end

