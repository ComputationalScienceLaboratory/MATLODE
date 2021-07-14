classdef (Abstract) Integrator < handle
    %INTEGRATOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Abstract, Constant)
        PartitionMethod
    end
    
    properties (SetAccess = immutable)
        Adaptive
        Datatype
    end
    
    methods (Abstract, Access = protected)
        [t, y, stats] = timeLoop(obj, f, tspan, y0, opts, varargin);
    end
    
    methods
        function obj = Integrator(adap, datatype)
            
            obj.Adaptive = adap;
            obj.Datatype = datatype;
        end
        
        function sol = integrate(obj, f, tspan, y0, varargin)
            
            [n,m] = size(tspan);
            
            if ~ismatrix(tspan)
                error('tspan cannot have this many dimensions');
            elseif n ~= 1 && m ~= 1
                error('tspan must be a vector')
            elseif (length(tspan) < 2)
                error('tspan must have a initial and final entry');
            elseif ~issorted(tspan, 'strictascend') && ~issorted(tspan, 'strictdescend')
                error('tspan must have unique values and be sorted in either descending or ascending order')
            end
            
            %Create options
            p = inputParser;
            opts = obj.matlodeSets(p, varargin{:});
            
            if isempty(opts.Dense) && length(tspan) > 2 && ~opts.StepSizeController.Adaptive
                error('IntegrateTo is not supported for a non-adaptable controller, please use DenseOutput instead')
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
                elseif isempty(f)
                    error('Partitioned system solves, require atleast one partition to be provided')
                else
                    %paritioning setup
                    if obj.PartitionMethod
                        [t, y, stats] = obj.timeLoop(f, tspan, y0, opts);
                    else
                       error('This method does not work with a partitioned f');
                    end
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
            p.addParameter('StepSizeController', matlode.stepsizecontroller.Fixed(1000));
            p.addParameter('ErrNorm', matlode.errnorm.InfNorm.errEstimate(sqrt(eps), sqrt(eps)));
            p.addParameter('Jacobian', []);
            p.addParameter('ChunkSize', 1000);
            p.addParameter('Dense', []);
            p.addParameter('MaxStep', inf);
            
            p.parse(varargin{:});
            
            opts = p.Results;
            
            %make sure a method can use a adaptive controller
            if ~obj.Adaptive && opts.StepSizeController.Adaptive
                error('Current method does not have the ability to use a adaptive step size controller');
            end
            
            
        end
        
    end
    
end

