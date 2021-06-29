classdef (Abstract) RungeKutta < matlode.OneStepIntegrator
    %RungeKutta integrator
    
    properties (Abstract, SetAccess = immutable)
        A
        B
        C
        BHat
        E
        Stage
        EmbeddedOrder
        Order
        FSAL
    end
    
    methods
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
    end
    
   
end

