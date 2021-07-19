classdef (Abstract) OneStepIntegrator < matlode.Integrator
    %One Step Integrator template
    
    
    methods (Access = protected)
        function obj = OneStepIntegrator(varargin)
            obj = obj@matlode.Integrator(varargin{:});
        end
        
        function opts = matlodeSets(obj, p, varargin)
            
            %OneStepIntegrator Specific options
            
            opts = matlodeSets@matlode.Integrator(obj, p, varargin{:});
            
        end
        
    end
end

