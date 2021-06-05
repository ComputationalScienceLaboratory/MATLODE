classdef (Abstract) Integrator < handle
    %INTEGRATOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Abstract, SetAccess = immutable)
        c
        adaptable
    end
    
    methods (Abstract)
        
        %integrate 
        sol = integrate(obj, func, tspan, y0, varargin)
        opts = matlodeSets(obj, varargin)
        
    end
    
    methods (Access = protected)
        
        
        function opts = matlodeSetsHelper(obj, p, varargin)
            
            p.addParameter('AbsTol', 1e-6);
            p.addParameter('RelTol', 1e-6);
            p.addParameter('InitialStep', 0);
            p.addParameter('StepSizeController', matlode.stepsizecontroller.Fixed(1000));
            p.addParameter('StepsN', 0);
            p.addParameter('ErrNorm', @(y, yHat, options) matlode.util.errnorm.standardNorm(y, yHat, options));
            p.addParameter('NewtonTol', 1e-7);
            p.addParameter('NewtonMaxItr', 25);
            p.addParameter('Jacobian', []);
            
            p.parse(varargin{:});
            
            opts = p.Results;
            
            %shorthanded fixed steps
            if opts.StepsN > 0
                opts.StepSizeController = matlode.stepsizecontroller.Fixed(opts.StepsN);
            end
            opts = rmfield(opts, 'StepsN');
            
            %make sure a method can use a adaptive controller
            if ~obj.adaptable && opts.StepSizeController.adaptive
                warning('Current method does not have the ability to use a adaptive step size controller');
                %overrider user choice???
            end
            
        end
    end
    
end

