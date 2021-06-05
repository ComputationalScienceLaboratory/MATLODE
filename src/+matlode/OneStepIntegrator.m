classdef (Abstract) OneStepIntegrator < matlode.Integrator
    %One Step Integrator template
    
    properties (Abstract, SetAccess = immutable)
        a
        b
        stage
    end
    
    methods (Access = protected)
        function opts = matlodeSetsHelper(obj, p, varargin)
            
            
            
            %OneStepIntegrator Specific options
            
            opts = matlodeSetsHelper@matlode.Integrator(obj, p, varargin{:});
            
        end
        
    end
end

