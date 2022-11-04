classdef (Abstract) Integrator < handle
    %INTEGRATOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Abstract, Constant)
        PartitionMethod
		PartitionNum
		MultirateMethod
    end
    
    properties (SetAccess = immutable)
        Adaptive
        Datatype
    end
    
    methods (Abstract, Access = protected)
        [t, y, stats] = timeLoop(obj, f, tspan, y0, opts);
		[t, y, stats] = timeLoopFixed(obj, f, tspan, y0, opts);
        [ynew, stages, stats] = timeStep(obj, f, t, y, dt, stages, prevAccept, stats);
        [err, stats] = timeStepErr(obj, f, t, y, ynew, dt, stages, ErrNorm, stats);
        [stages, stats] = timeLoopBeforeLoop(obj, f, f0, t0, y0, stats);
        [q] = timeLoopInit(obj);
    end
    
    methods
		%Integrate over a adaptive step
        function obj = Integrator(adap, datatype)
            
            obj.Adaptive = adap;
            obj.Datatype = datatype;
        end
        
        function sol = integrate(obj, f, tspan, y0, varargin)

			if ~obj.Adaptive
				error('Method does not support adaptive step control. Call integrateFixed for integration')
			end
            
            [n,m] = size(tspan);
            
            if (n == 1 && m == 1)
                error('tspan must be a vector');
			elseif (n > 1 && m > 1)
                error('tspan must be a vector');
            elseif ~issorted(tspan, 'strictascend') && ~issorted(tspan, 'strictdescend')
                error('tspan must have unique values and be sorted in either descending or ascending order')
            end
            
            %Create options

			p = inputParser();
			p.KeepUnmatched = true;
            opts = obj.matlodeSets(p, varargin{:});
            

			fmod = matlode.Model(f, varargin{:});

            [t, y, stats] = obj.timeLoop(fmod, tspan, y0, opts);

            
            sol.t = t;
            sol.y = y;
            sol.stats = stats;

			sol.model = fmod;
		end

		%Integrate over a fixed interval
		function sol = integrateFixed(obj, f, tspan, y0, varargin)
            
            [n,m] = size(tspan);
            
            if (n == 1 && m == 1)
                error('tspan must be a vector');
			elseif (n > 1 && m > 1)
                error('tspan must be a vector');
            elseif ~issorted(tspan, 'strictascend') && ~issorted(tspan, 'strictdescend')
                error('tspan must have unique values and be sorted in either descending or ascending order')
            end
            
            %Create options
			p = inputParser();
			p.KeepUnmatched = true;
            opts = obj.matlodeSets(p, varargin{:});

			fmod = matlode.Model(f, varargin{:});

            [t, y, stats] = obj.timeLoop(fmod, tspan, y0, opts);
            
            sol.t = t;
            sol.y = y;
            sol.stats = stats;

			sol.model = fmod;
		end
    end
    
    methods (Access = protected)
        
        function opts = matlodeSets(obj, p, varargin)
            
            %Assign to controllers
			p.addParameter('OTPPath', []);
            p.addParameter('StepSizeController', matlode.stepsizecontroller.StandardController);
            p.addParameter('ErrNorm', matlode.errnorm.InfNorm(sqrt(eps), sqrt(eps)));
            p.addParameter('ChunkSize', 1000);
			p.addParameter('FullTrajectory', false);
            p.addParameter('Dense', []);
            p.addParameter('MaxStep', inf);
			
            
            p.parse(varargin{:});
            
            opts = p.Results;
            
            
        end
        
    end
    
end

