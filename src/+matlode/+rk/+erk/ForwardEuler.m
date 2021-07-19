classdef ForwardEuler < matlode.rk.erk.ERK
    methods
        function obj = ForwardEuler(datatype)
            %TODO
            %include reference
            
            %p = 1 & hatp = 0 & s = 1
            
            if nargin == 0
                datatype = 'double';
            end
            
            caster = @(x) matlode.util.CoefficentTransformers.transform(x,datatype);
            
            a = caster('[[0]]');
             
            b = caster('[1]');
            
            bHat = [];
            
            e = [];
            
            c = caster('[0]');
            
            order = 1;
            
            embbededOrder = 0;
            
            obj = obj@matlode.rk.erk.ERK(a, b, bHat, c, e, order, embbededOrder);
            
            %p* = 1 & s* = 1
        end
        
        function sol = integrate(obj, f, tspan, y0, varargin)
            sol = integrate@matlode.rk.erk.ERK(obj, f, tspan, y0,'Dense', matlode.denseoutput.Linear(obj.B), 'StepSizeController', matlode.stepsizecontroller.Fixed(1000), varargin{:});
        end
    end
end

