classdef PID < matlode.stepsizecontroller.SoderlandController
    % Gustaffson "Control-Theoretic Techniques for Stepsize Selection in
    % Implicit Runge-Kutta Methods"
    % Solving ODEs II Chapter IV pg 124
    methods
        function obj = PID(varargin)
            
            a = [-2, 1];
            b = [1];
            
            AQfunc = @(q, A) (A / q);
            
            obj = obj@matlode.stepsizecontroller.SoderlandController('A', a, 'B', b, 'AQfunc', AQfunc, varargin{:});
        end
    end
end

