classdef RODAS3 < matlode.rosenbrock.Rosenbrock
    %RODAS3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = immutable)
        Dense
    end
    
    methods
        function obj = RODAS3(datatype)
            if nargin == 0
                datatype = 'double';
            end
            
            
            caster = @(x) matlode.util.CoefficentTransformers.transform(x,datatype);
			
            gammadia = caster('[1/2, 1/2, 1/2, 1/2]');
            alphasum = caster('[0, 0, 1, 1]');
            gammasum = caster('[1/2, 3/2, 0, 0]');
            
            a = caster(join(['[[0, 0, 0, 0];',...
                             '[0, 0, 0, 0];',...
                             '[2, 0, 0, 0];',...
                             '[2, 0, 1, 0]]']));
            
            c = caster(join(['[[0, 0, 0, 0];',...
                             '[4, 0, 0, 0];',...
                             '[1, -1, 0, 0];',...
                             '[1, -1, -8/3, 0]]']));
            m = caster('[2, 0, 1, 1]');
            me = caster('[0, 0, 0, 1]');
            
            order = 3;
            
            embbededOrder = 2;
            
            obj = obj@matlode.rosenbrock.Rosenbrock(gammadia, gammasum, alphasum, a, c, m, me, order, embbededOrder);
            
            
            obj.Dense = matlode.denseoutput.Linear(m);
        end
        
        function sol = integrate(obj, f, tspan, y0, varargin)
            sol = integrate@matlode.rosenbrock.Rosenbrock(obj, f, tspan, y0, 'StepSizeController', matlode.stepsizecontroller.StandardController, 'Dense', obj.Dense, varargin{:});
        end
    end
end

