classdef BackwardEuler < matlode.rk.dirk.DIRK
	%BackwardEuler
	% p = 1, s = 1

	methods
		function obj = BackwardEuler(datatype)
			arguments
				datatype(1,1) string = 'double';
			end

			caster = @(x) matlode.util.CoefficentTransformers.transform(x,datatype);


			a = caster('[[1]]');

			b =  caster('[1]');

			bHat = [];

			e = [];

			c = caster('[1]');

			order = 1;

			embbededOrder = 0;

			obj = obj@matlode.rk.dirk.DIRK(a, b, bHat, c, e, order, embbededOrder);

		end

		function sol = integrate(obj, f, tspan, y0, varargin)
			sol = integrate@matlode.rk.dirk.DIRK(obj, f, tspan, y0, 'StepSizeController', matlode.stepsizecontroller.StandardController, 'Dense', obj.DenseOut, varargin{:});
		end
	end
end

