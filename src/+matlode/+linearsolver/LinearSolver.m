classdef (Abstract) LinearSolver
    methods (Abstract)
        system = preprocess(obj, updateState, updateTimestep, oldSystem, t, y, f, m1, mass, m2, jac);
        
        sol = solve(obj, x, system, t, y, f, m1, mass, m2, jac);
    end
end

