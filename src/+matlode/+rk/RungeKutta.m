classdef (Abstract) RungeKutta < matlode.OneStepIntegrator
    %RungeKutta integrator
    
    properties (SetAccess = immutable)
        A
        B
        C
        BHat
        E
        StageNum
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
            obj.StageNum = size(b, 2);
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
       
        
        function [stages, stats] = timeLoopBeforeLoop(obj, f, f0, t0, y0, stats)
            
            stages = zeros(length(y0), obj.StageNum);
            
            if obj.FSAL
                if isempty(f0)
                    stages(:, end) = f.F(t0, y0);
                    stats.FEvals = stats.FEvals + 1;
                else
                    stages(:, end) = f0;
                end
            end
        end
        
        function [q] = timeLoopInit(obj)
            q = min(obj.Order, obj.EmbeddedOrder);
        end
    end
    
   
end

