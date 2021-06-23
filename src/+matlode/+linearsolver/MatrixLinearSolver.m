classdef MatrixLinearSolver < matlode.linearsolver.LinearSolver
    properties (SetAccess = immutable, GetAccess = private)
        Solver
    end
    
    methods
        function obj = MatrixLinearSolver(solver)
            if nargin < 1
                solver = @mldivide;
            end
            
            obj.Solver = solver;
        end
        
        function system = preprocess(~, updateState, updateTimestep, system, t, y, ~, m1, mass, m2, jac)
            init = isempty(system);
            
            if updateState
                if ~isnumeric(jac)
                    system.jac = jac(t, y);
                elseif init
                    system.jac = jac;
                end
                
                if ~isnumeric(mass)
                    system.mass = mass(t, y);
                elseif init
                    if isempty(mass)
                        system.mass = eye(size(system.jac), 'like', system.jac);
                    else
                        system.mass = mass;
                    end
                end
            end
            
            if updateState || updateTimestep
                system.matrix = kron(m1, system.mass);
                if m2 ~= 0
                    system.matrix = system.matrix - kron(m2, jac);
                end
            end
        end
        
        function sol = solve(obj, x, system, ~, ~, ~, ~, ~, ~, ~)
            sol = obj.Solver(system.matrix, x);
        end
    end
end

