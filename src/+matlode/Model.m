classdef Model < handle
	%Model describes the model to propagate forward in time
	%Will contain all properties that are associated with the model such as
	%F,Jacobian, Mass Matrix, etc..
	%All time Integration methods are passed this object
	
	properties
		F
		AddPartitionNum
		ComponentPartitionNum

		Jacobian
        JacobianVectorProduct
		JacobianAdjointVectorProduct
        JPattern
        Vectorized

        Mass
		MVectorProduct

		%For complatency with MATLAB standard. May have DAEModel instead
        MassSingular
        MStateDependence
        MvPattern

		PartialDerivativeTime

		PartialDerivativeParameters
		
		%FIXME: Not Supported
		Events
		OnEvent

        NonNegative
	end

	properties(GetAccess = protected, Constant)
		%% Constant properties for conversions
		MATLABSetVars = {'Jacobian', 'JPattern', 'Vectorized', 'Mass', 'MassSingular', 'MStateDependence', ...
			'MvPattern', 'MassSingular', 'Events', 'NonNegative'}
		OTPSetVars = {'Jacobian', 'JPattern', 'Vectorized', 'Mass', 'MStateDependence', ...
			'MvPattern', 'MassSingular', 'Events', 'F', 'JacobianVectorProduct', 'JacobianAdjointVectorProduct', ...
			'PartialDerivativeTime' ,'PartialDerivativeParameters', ...
			'OnEvent', 'NonNegative'}
	end
	
	methods
		function obj = Model(f,varargin)
	
			%Handle varargin
    		p = inputParser;
    		p.KeepUnmatched = true;

			%Handle what f could be
			if isa(f, 'function_handle') || (isa(f, 'cell') && cellfun(@(x) isa(x, 'function_handle'), f))
				%% Create new model based of f
				obj.F = f;
	
				p.addParameter('Jacobian', []);
				p.addParameter('JPattern', []);
				p.addParameter('JacobianVectorProduct', []);
				p.addParameter('Vectorized', []);
	
				p.addParameter('Mass', []);
				p.addParameter('MassSingular', false);
				p.addParameter('MStateDependence', false);
				p.addParameter('MvPattern', []);
	
				p.addParameter('PartialDerivativeTime', []);
	
				p.addParameter('PartialDerivativeParameters', []);
	
	
				p.addParameter('Events', []);
				p.addParameter('OnEvent', []);
	
				p.addParameter('NonNegative', []);
			elseif isa(f, 'otp.RHS')
				%% OTP RHS Parsing
				obj = otpRHSToModel(obj, f);

				for i = 1:length(obj.OTPSetVars)
					p.addParameter(obj.OTPSetVars{i}, obj.(obj.OTPSetVars{i}));
				end

			elseif isa(f, 'matlode.model.Model')
				%% Model Copying
				
				obj = modelCopy(obj, f);

				props = properties(obj);

				for i = 1:length(props)
					p.addParameter(props{i}, obj.(props{i}));
				end
			end

    		%Handles the varargin
			%Will overwrite model with new parameters
    		p.parse(varargin{:});

			parms = p.Results;
	
			fields = fieldnames(parms);
			props = properties(obj);
			for i = 1:length(fields)
				f = fields{i};
				index = find(strcmp(props, f), 1);
				if ~isempty(index)
					obj.(f) = parms.(f);
				end
			end

			%Check Dimension of F
			%TODO Turn into Function
			if isa(f, 'cell')
				obj.AddPartitionNum = size(f,2);
				obj.ComponentPartitionNum = size(f,1);
			else
				obj.AddPartitionNum = 1;
				obj.ComponentPartitionNum = 1;
			end

		end

		function obj = setF(obj, f)
			obj.F = f;
		end
		
		%Take a preexsisting model and copy into a new one
		function obj = modelCopy(obj, otherModel)

			props = properties(obj);

			for i = 1:length(props)
				obj.(props{i}) = otherModel.(props{i});
			end
		end

		%Convert Model to matlab ode set
		function [f, odeopts] = modelToMATLABOdeSets(obj)
			f = obj.F;
			odeopts = odeset([]);
			
			for i = 1:length(obj.MATLABSetVars)
				odeopts.(obj.MATLABSetVars{i}) = obj.(obj.MATLABSetVars{i});
			end

		end

		%Convert matlab ode set to Model
		function obj = matlabOdeSetsToModel(obj, f, odeopts)
			obj.F = f;

			for i = 1:length(obj.MATLABSetVars)
				obj.(obj.MATLABSetVars{i}) = odeopts.(obj.MATLABSetVars{i});
			end
		end

		%Convert  OTP RHS object to Model
		function obj = otpRHSToModel(obj, rhs)
			
			for i = 1:length(obj.OTPSetVars)
				obj.(obj.OTPSetVars{i}) = rhs.(obj.OTPSetVars{i});
			end
		end

	end

end

