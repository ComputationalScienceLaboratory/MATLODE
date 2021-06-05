classdef (Abstract) RungeKutta < matlode.OneStepIntegrator
    %RungeKutta integrator
    
    properties (Abstract, SetAccess = immutable)
        bHat
        embeddedOrder
        order
        FSAL
    end
    
    methods
        function obj = fromCoeffs(a, b, bHat, c, order, embededOrder)
            
            
            if istril(a)
                if any(diag(a))
                    %TODO
                    % set as DIRK
                    return
                end
                
                obj = ERK(a, b, bHat, c, order, embededOrder);
                return
            end
            
            %TODO
            
            %Set as FIRK
                
        end
    end
    
    methods (Access = protected)
        function opts = matlodeSetsHelper(obj, p, varargin)
            
            %RungeKutta sepcific options
            
            opts = matlodeSetsHelper@matlode.OneStepIntegrator(obj, p, varargin{:});
            
        end
    end
    
   
end

