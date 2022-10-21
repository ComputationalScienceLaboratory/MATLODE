classdef LinBackEuler < matlode.rosenbrock.Rosenbrock
    %LINBACKEULER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Dense
    end
    
    methods
        function obj = LinBackEuler(varargin)
            
            p = inputParser;
            p.KeepUnmatched = true;
            
            p.addParameter('datatype', 'double');
            
            p.parse(varargin{:});
            opts = p.Results;

            datatype = opts.datatype;
            
            caster = @(x) matlode.util.CoefficentTransformers.transform(x,datatype);
            
            alpha = caster('[0]');
            gamma = caster('[1]');
            b = caster('[1]');
            be = [];

			[gammadia, gammasum, alphasum, a, c, m, me] = RosCoefMethTrans(gamma, alpha, b, be);
            
            order = 1;
            
            embbededOrder = 0;
            
            
            obj = obj@matlode.rosenbrock.Rosenbrock(gammadia, gammasum, alphasum, a, c, m, me, order, embbededOrder);
            
            obj.Dense = matlode.denseoutput.Linear(b);
            
        end
        
        function sol = integrate(obj, f, tspan, y0, varargin)
            sol = integrate@matlode.rosenbrock.Rosenbrock(obj, f, tspan, y0, 'StepSizeController', matlode.stepsizecontroller.Fixed(1000), 'Dense', obj.DenseOut, varargin{:});
        end
    end
end

