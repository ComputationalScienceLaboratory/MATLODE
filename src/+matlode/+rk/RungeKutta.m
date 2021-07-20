classdef (Abstract) RungeKutta < matlode.OneStepIntegrator
    %RungeKutta integrator
    
    properties (SetAccess = immutable)
        A
        B
        C
        BHat
        E
        Stage
        EmbeddedOrder
        Order
        FSAL
        FsalStart
    end
    
    methods
        function obj = RungeKutta(a, b, bHat, c, e, order, embeddedOrder)
            
            obj = obj@matlode.OneStepIntegrator(~isempty(bHat), class(a));
            
            obj.A = a;
            obj.B = b;
            obj.BHat = bHat;
            obj.C = c;
            obj.E = e;
            obj.EmbeddedOrder = embeddedOrder;
            obj.Order = order;
            obj.Stage = size(b, 2);
            obj.FSAL = all(a(end, :) == b) && all(a(1, :) == 0);
            obj.FsalStart = uint32(obj.FSAL) + 1;
        end
        
        function obj = fromCoeffs(a, b, bHat, c, e, bTheta, order, embededOrder)
            
            
            if istril(a)
                if any(diag(a))
                    %TODO
                    % set as DIRK
                    return
                end
                
                obj = ERK(a, b, bHat, c, e, bTheta, order, embededOrder);
                return
            end
            
            %TODO
            
            %Set as FIRK
                
        end
    end
    
    methods (Access = protected)
        function opts = matlodeSets(obj, p, varargin)
            
            %RungeKutta sepcific options
            
            opts = matlodeSets@matlode.OneStepIntegrator(obj, p, varargin{:});
            
        end
       
        
        function [k, fevals, fevalIterCounts] = timeLoopBeforeLoop(obj, f0, t0, y0)
            
            k = zeros(length(y0), obj.Stage);
            
            fevalIterCounts = double(obj.Stage - obj.FsalStart + 1);
            
            if obj.FSAL
                if isempty(f0)
                    k(:, end) = f(t0, y0);
                    fevals = 1;
                else
                    k(:, end) = f0;
                    fevals = 0;
                end
            end
        end
        
        function [q] = timeLoopInit(obj)
            q = min(obj.Order, obj.EmbeddedOrder);
        end
    end
    
   
end

