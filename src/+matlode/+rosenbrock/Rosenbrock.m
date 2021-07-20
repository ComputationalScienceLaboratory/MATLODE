classdef Rosenbrock < matlode.OneStepIntegrator
    %ROSENBROCK Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = immutable)
        A
        B
        C
        G
        BHat
        E
        Stage
        EmbeddedOrder
        Order
        FSAL
        FsalStart
    end
    
    methods
        function obj = Rosenbrock(a, b, bHat, c, e, g, order, embeddedOrder)
            obj.A = a;
            obj.B = b;
            obj.BHat = bHat;
            obj.C = c;
            obj.E = e;
            obj.G = g;
            obj.EmbeddedOrder = embeddedOrder;
            obj.Order = order;
            obj.Stage = size(b, 2);
            obj.FSAL = all(a(end, :) == b) && all(a(1, :) == 0);
            obj.FsalStart = uint32(obj.FSAL) + 1;
        end
        
        
    end
    
    methods (Access = protected)
        function [k, yN, err] = timeStep(obj,  f, t, y, h, k, ErrNorm, prevAccept)
            
        end
        
        function [k, yN, err] = timeStepErr(obj,  f, t, y, h, k, ErrNorm, prevAccept)
            
        end
        
        function [k, fevals, fevalIterCounts] = timeLoopBeforeLoop(obj, f0, t0, y0)
            
        end
        
        function [q] = timeLoopInit(obj)
            
        end
    end
end

