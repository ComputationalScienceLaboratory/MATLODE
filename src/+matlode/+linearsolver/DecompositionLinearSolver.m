classdef DecompositionLinearSolver < matlode.linearsolver.MatrixLinearSolver
    properties (SetAccess = immutable, GetAccess = private)
        DecompositionType
		DecompositionArgs
    end
    
    methods
		function obj = DecompositionLinearSolver(type, solver, decompArgs)
			arguments
				type(1,1) string = 'lu';
				solver(1,1) function_handle = @mldivide;
				decompArgs(1,:) cell = {};
			end
            obj = obj@matlode.linearsolver.MatrixLinearSolver(solver, {});
			obj.DecompositionType = type;
			obj.DecompositionArgs = decompArgs;
        end
        
        function [stats] = preprocess(obj, f, t, y, reeval, mass_scale, jac_scale, stats)
            
            [stats] = preprocess@matlode.linearsolver.MatrixLinearSolver(obj, f, t, y, reeval, mass_scale, jac_scale, stats);
            
			if isempty(obj.DecompositionArgs)
				obj.system = decomposition(obj.system, obj.DecompositionType);
			else
				obj.system = decomposition(obj.system, obj.DecompositionType, obj.DecompositionArgs);
			end
            
            stats.nDecompositions = stats.nDecompositions + 1;
        end
    end
end

