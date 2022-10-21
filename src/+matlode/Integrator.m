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
            
            t = [];
            y = [];
            stats = [];

			new_model = false;
            
            %Compute
            if isa(f, 'matlode.Model')
				if isempty(f.Mass)
					f.Mass = eye(length(y0));
				end
                
                [t, y, stats] = obj.timeLoop(f, tspan, y0, opts);
			elseif isa(f, 'function_handle') || isa(f, 'cell')
				if isa(f, 'cell') && ~obj.PartitionMethod
					error('This method does not support partitioning')
				end
				f_mod = matlode.Model(f, varargin{:});
				new_model = true;
				if isempty(f_mod.Mass)
					f_mod.Mass = eye(length(y0));
				end

                [t, y, stats] = obj.timeLoop(f_mod, tspan, y0, opts);

			elseif isa(f, 'otp.RHS')
				f_mod = matlode.Model(f, varargin{:});
				new_model = true;
				if isempty(f_mod.Mass)
					f_mod.Mass = eye(length(y0));
				end

                [t, y, stats] = obj.timeLoop(f_mod, tspan, y0, opts);
                
            else
                error('f must be a ''function_handle'', ''cell'', OTP RHS, or Model object');
            end
            
            sol.t = t;
            sol.y = y;
            sol.stats = stats;

			if new_model
				sol.model = f_mod;
			else
				sol.model = f;
			end
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
            
            
            t = [];
            y = [];
            stats = [];

			new_model = false;
            
            %Compute
            if isa(f, 'matlode.Model')
                
                [t, y, stats] = obj.timeLoopFixed(f, tspan, y0, opts);
			elseif isa(f, 'function_handle') || isa(f, 'cell')
				if isa(f, 'cell') && ~obj.PartitionMethod
					error('This method does not support partitioning')
				end
				f_mod = matlode.Model(f, varargin{:});
				new_model = true;

                [t, y, stats] = obj.timeLoopFixed(f_mod, tspan, y0, opts);

			elseif isa(f, 'otp.RHS')
				f_mod = matlode.Model(f, varargin{:});
				new_model = true;

                [t, y, stats] = obj.timeLoopFixed(f_mod, tspan, y0, opts);
                
            else
                error('f must be a ''function_handle'', ''cell'', OTP RHS, or Model object');
            end
            
            sol.t = t;
            sol.y = y;
            sol.stats = stats;

			if new_model
				sol.model = f_mod;
			else
				sol.model = f;
			end
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

