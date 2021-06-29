classdef (Abstract) OneStepIntegrator < matlode.Integrator
    %One Step Integrator template
    
    
    methods (Access = protected)
        function opts = matlodeSets(obj, p, varargin)
            
            %OneStepIntegrator Specific options
            
            opts = matlodeSets@matlode.Integrator(obj, p, varargin{:});
            
        end
        
    end
end

