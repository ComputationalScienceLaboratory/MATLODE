classdef PI < matlode.stepsizecontroller.SoderlandController
    %PI controller
    % From Solving ODE II Chapter IV.2
    
    methods
        function obj = PI(varargin)
            a = [-0.7, 0.4];
            
            AQfunc = @(q, A) (A / q);
            
            obj = obj@matlode.stepsizecontroller.SoderlandController('A', a, 'AQfunc', AQfunc, varargin{:});
        end
    end
end

