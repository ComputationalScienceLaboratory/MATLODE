classdef DecompositionLinearSolver < matlode.linearsolver.MatrixLinearSolver
    properties (SetAccess = immutable, GetAccess = private)
        Args
    end
    
    methods
        function obj = DecompositionLinearSolver(varargin)
            obj = obj@matlode.linearsolver.MatrixLinearSolver();
            obj.Args = varargin;
        end
        
        function system = preprocess(obj, updateState, updateTimestep, system, t, y, f, m1, mass, m2, jac)
            system = preprocess@matlode.linearsolver.MatrixLinearSolver( ...
                obj, updateState, updateTimestep, system, t, y, f, m1, mass, m2, jac);
            
            if updateState || updateTimestep
                system.matrix = decomposition(system.matrix, obj.Args{:});
            end
        end
    end
end

