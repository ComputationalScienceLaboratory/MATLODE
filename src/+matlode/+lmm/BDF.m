classdef BDF < matlode.Integrator
	methods
		function obj = BDF(datatype)
			arguments
				datatype(1,1) string = 'double';
			end

			obj = obj@matlode.Integrator(true, datatype)
		end
	end

	methods (Access = protected)
		function opts = matlodeSets(obj, p, varargin)

			%BDF Specific options

			opts = matlodeSets@matlode.Integrator(obj, p, varargin{:});

		end
	end
end