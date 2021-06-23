classdef MatrixFreeLinearSolver < matlode.linearsolver.LinearSolver
    properties (SetAccess = immutable, GetAccess = private)
        Solver
    end
    
    methods
        function obj = MatrixFreeLinearSolver(solver)
            if nargin < 1
                solver = @gmres;
            end
            
            obj.Solver = solver;
        end
        
        function system = preprocess(~, ~, ~, system, ~, ~, ~, ~, ~, ~, ~)
        end
        
        function sol = solve(obj, x, ~, t, y, ~, m1, mass, m2, jac)
            if length(m1) == 1
                if isempty(jac)
                    jvp = @(v) m1 * mass(t, y, v);
                elseif isempty(mass)
                    jvp = @(v) m1 * v - m2 * jac(t, y, v);
                else
                    jvp = @(v) m1 * mass(t, y, v) - m2 * jac(t, y, v);
                end
            else
                error('Matrix-free block systems are not supported yet');
            end
            sol = obj.Solver(jvp, x);
        end
    end
end

