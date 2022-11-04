classdef Model < handle
	%Model describes the model to propagate forward in time
	%Will contain all properties that are associated with the model such as
	%F,Jacobian, Mass Matrix, etc..
	%All time Integration methods are passed this object
	
	properties
		F
		FOperatorType

		Jacobian
		JLinearOperatorType
        JacobianVectorProduct
		JacobianAdjointVectorProduct
        JPattern
        Vectorized

        Mass
		MLinearOperatorType
		MVectorProduct
        MassSingular
        MStateDependence
        MvPattern

		PartialDerivativeTime
		PDTOperatorType

		PartialDerivativeParameters
		PDPOperatorType

        HessianVectorProduct
        HessianAdjointVectorProduct
		
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
			'PartialDerivativeTime' ,'PartialDerivativeParameters', 'HessianVectorProduct', 'HessianAdjointVectorProduct', ...
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
				
				p.addParameter('FOperatorType', matlode.OperatorType.TSDep);
	
				p.addParameter('Jacobian', []);
				p.addParameter('JLinearOperatorType', matlode.LinearOperatorType.Empty);
				p.addParameter('JPattern', []);
				p.addParameter('JacobianVectorProduct', []);
				p.addParameter('Vectorized', []);
	
				p.addParameter('Mass', []);
				p.addParameter('MLinearOperatorType', matlode.LinearOperatorType.Empty);
				p.addParameter('MassSingular', false);
				p.addParameter('MStateDependence', false);
				p.addParameter('MvPattern', []);
	
				p.addParameter('PartialDerivativeTime', []);
				p.addParameter('PDTOperatorType', matlode.OperatorType.Zero);
	
				p.addParameter('PartialDerivativeParameters', []);
				p.addParameter('PDPOperatorType', matlode.OperatorType.Zero);
	
				p.addParameter('HessianVectorProduct', []);
				p.addParameter('HessianAdjointVectorProduct', []);
	
	
				p.addParameter('Events', []);
				p.addParameter('OnEvent', []);
	
				p.addParameter('NonNegative', []);
			elseif isa(f, 'otp.RHS')
				%% OTP RHS Parsing
				obj = otpRHSToModel(obj, f);

				for i = 1:length(obj.OTPSetVars)
					p.addParameter(obj.OTPSetVars{i}, obj.(obj.OTPSetVars{i}));
				end
        		
				p.addParameter('PDTOperatorType', matlode.OperatorType.Zero);
				p.addParameter('PDPOperatorType', matlode.OperatorType.Zero);
				p.addParameter('MLinearOperatorType', matlode.LinearOperatorType.Empty);
				p.addParameter('JLinearOperatorType', matlode.LinearOperatorType.Empty);
				p.addParameter('FOperatorType', matlode.OperatorType.TSDep);

			elseif isa(f, 'matlode.Model')
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

			obj = modelOperatorTypeSet(obj);

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

	methods (Access = protected)

		function obj = modelOperatorTypeSet(obj)
			if isnumeric(obj.Jacobian)
				obj.JLinearOperatorType = matlode.LinearOperatorType.Constant;
			end

			%% Operator Type Checking
			%TODO: Support Constant Operators

			obj = obj.opTypeSet('PartialDerivativeTime', 'PDTOperatorType', matlode.OperatorType.Zero, matlode.OperatorType.TSDep);
			obj = obj.opTypeSet('Jacobian', 'JLinearOperatorType', matlode.LinearOperatorType.Empty, matlode.LinearOperatorType.TSDepedent);
			if obj.JLinearOperatorType == matlode.LinearOperatorType.Empty && ~isempty(obj.JacobianVectorProduct)
				obj.JLinearOperatorType = matlode.LinearOperatorType.TSDepedent;
			end
			% Currently only support constant Mass
			% TODO: Support Time Dependent Mass
			% TODO: Support State Dependent Mass
			obj = obj.opTypeSet('Mass', 'MLinearOperatorType', matlode.LinearOperatorType.Empty, matlode.LinearOperatorType.Constant);
			if obj.MLinearOperatorType == matlode.LinearOperatorType.Empty && ~isempty(obj.MVectorProduct)
				obj.MLinearOperatorType = matlode.LinearOperatorType.TSDepedent;
			end

			switch(obj.MLinearOperatorType)
				case matlode.LinearOperatorType.Zero
				case matlode.LinearOperatorType.Identity
				case matlode.LinearOperatorType.Constant
				otherwise
					if strcmp(obj.MassSingular, 'yes')
						error('Non-Constant Singular Mass Matrices are not Supported')
					end
			end
		end

		function obj = opTypeSet(obj, op, optype, opttypedef, optnew)
			if ~isempty(obj.(op)) && (obj.(optype) == opttypedef)
				obj.(optype) = optnew;
			end

		end
	end
end

