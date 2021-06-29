classdef ClassicRK4 < matlode.rk.erk.ERK
    properties (SetAccess = immutable)
        DenseOut
    end
    
    methods
        function obj = ClassicRK4(datatype)
            
            %TODO
            %include reference
            
            % p = 4 & hatp = 0 & s = 4
            
            if nargin == 0
                datatype = 'double';
            end
            
            caster = @(x) matlode.util.CoefficentTransformers.transform(x,datatype);
            
            a = caster(join(['[[0,   0,   0, 0];',...
                             '[1/2, 0,   0, 0];',...
                             '[0,   1/2, 0, 0];',...
                             '[0,   0,   0, 1]]']));
             
            b = caster('[1/6, 1/3, 1/3, 1/6]');
            
            bHat = [];
            
            e = [];
            
            c = caster('[0, 1/2, 1/2, 1]');
            
            order = 4;
            
            embbededOrder = 0;
            
            obj = obj@matlode.rk.erk.ERK(a, b, bHat, c, e, order, embbededOrder);
            
            %p* = 3 & s* = 3
            denseMatrix = caster(join(['[1, -3/2, 2/3];',...
                                       '[0, 1,    -2/3]',...
                                       '[0, 1,    -2/3]',...
                                       '[0, 1/2,  2/3]']));
            
            obj.DenseOut = matlode.denseoutput.RKDense(denseMatrix);
        end
        
        function sol = integrate(obj, f, tspan, y0, varargin)
            sol = integrate@matlode.rk.erk.ERK(obj, f, tspan, y0, 'StepSizeController', matlode.stepsizecontroller.Fixed(1000), varargin{:});
        end
    end
end

