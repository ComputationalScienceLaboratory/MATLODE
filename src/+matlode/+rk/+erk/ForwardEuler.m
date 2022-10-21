classdef ForwardEuler < matlode.rk.erk.ERK
%Method: Forward Euler
% p = 1 s = 1 pe = -1
%Reference: 

    methods
        function obj = ForwardEuler(datatype)
            arguments
				datatype(1,1) string = 'double';
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
            sol = integrate@matlode.rk.erk.ERK(obj, f, tspan, y0,'Dense', obj.DenseOut, varargin{:});
        end
    end
end

