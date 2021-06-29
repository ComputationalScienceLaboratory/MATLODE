classdef (Abstract) Integrator < handle
    %INTEGRATOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Abstract, Constant)
        PartitionMethod
    end
    
    properties (Abstract, SetAccess = immutable)
        Adaptive
        Datatype
    end
    
    methods (Abstract, Access = protected)
        [t, y, stats] = timeLoop(obj, f, tspan, y0, opts, varargin);
    end
    
    methods
        function sol = integrate(obj, f, tspan, y0, varargin)
            
            if(length(tspan) < 2)
                error('Unsupported tspan entry');
            end
            
            if ~issorted(tspan) && ~issorted(tspan, 'descend')
                error('tspan must be sorted in order')
            end
            
            %TODO
            %optimize options for when only options struct is given
            
            %Create options
            p = inputParser;
            opts = obj.matlodeSets(p, varargin{:});
            
            if isempty(opts.Dense) && length(tspan) > 2 && isa(opts.StepSizeController, 'matlode.stepsizecontroller.Fixed')
                error('IntegrateTo is not supported for fixed step size, please use DenseOutput instead')
            end
            
            t = [];
            y = [];
            stats = [];
            
            %Compute
            if isa(f, 'function_handle') && ~obj.PartitionMethod
                
                [t, y, stats] = obj.timeLoop(f, tspan, y0, opts);
                
            elseif isa(f, 'cell') && cellfun(@(x) isa(x, 'function_handle'), f)
                
                %determine size of partition
                if length(f) == 1 && ~obj.PartitionMethod
                    [t, y, stats] = obj.timeLoop(f{1}, tspan, y0, opts);
                elseif length(f) > 1
                    %TODO
                    %paritioning setup
                    
                else
                    error('Need at least one function provided')
                end
                
            else
                error('f must be a ''function_handle'' or a ''cell''');
            end
            
            sol.t = t;
            sol.y = y;
            sol.stats = stats;
        end
    end
    
    methods (Access = protected)
        
        
        function opts = matlodeSets(obj, p, varargin)
            
            %Assign to controllers
            p.addParameter('InitialStep', 0);
            p.addParameter('StepSizeController', matlode.stepsizecontroller.Fixed(1000));
            p.addParameter('ErrNorm', matlode.errnorm.HistNorm.errEstimate(sqrt(eps), sqrt(eps)));
            p.addParameter('Jacobian', []);
            p.addParameter('ChunkSize', 1000);
            p.addParameter('Dense', [] );
            p.addParameter('MaxStep', inf);
            
            p.parse(varargin{:});
            
            opts = p.Results;
            
            %make sure a method can use a adaptive controller
            if ~obj.Adaptive && opts.StepSizeController.Adaptive
                warning('Current method does not have the ability to use a adaptive step size controller');
                %overrider user choice???
            end
            
            
        end
        
    end
    
end

