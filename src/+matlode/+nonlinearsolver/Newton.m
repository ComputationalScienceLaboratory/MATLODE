classdef Newton < matlode.nonlinearsolver.NonlinearSolver
	%NEWTON Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
		Property1
	end
	
	methods
		function obj = Newton(linsolve, args)
			arguments
				linsolve(1,1) matlode.linearsolver.LinearSolver = {};
				args(1,:) cell = {};
			end
            
            obj = obj@matlode.nonlinearsolver.NonlinearSolver(linsolve, args{:});
		end
		
		function outputArg = method1(obj,inputArg)
			%METHOD1 Summary of this method goes here
			%   Detailed explanation goes here
			outputArg = obj.Property1 + inputArg;
		end
	end
end

