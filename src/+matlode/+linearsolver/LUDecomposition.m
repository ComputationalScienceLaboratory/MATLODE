classdef LUDecomposition < matlode.linearsolver.MatrixLinearSolver
    
    methods
		function obj = LUDecomposition()
            obj = obj@matlode.linearsolver.MatrixLinearSolver(@mldivide, {});
        end
        
        function [stats] = preprocess(obj, f, t, y, reeval, mass_scale, jac_scale, stats)
            
            [stats] = preprocess@matlode.linearsolver.MatrixLinearSolver(obj, f, t, y, reeval, mass_scale, jac_scale, stats);
            
			if issparse(obj.system)
				[L,U,P,Q,D] = lu(obj.system);
				obj.system = @(x) Q*(U\(L\(P*(D \ x))));
			else
				[L,U,P] = lu(obj.system);
				obj.system = @(x) U\(L\(P*x));
			end
            
            stats.nDecompositions = stats.nDecompositions + 1;
		end

		function [sol, stats] = solve(obj, x, stats)
            sol = obj.system(x);
            
            stats.nLinearSolves = stats.nLinearSolves + 1;
        end
    end
end

