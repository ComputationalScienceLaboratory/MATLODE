classdef H321 < matlode.stepsizecontroller.SoderlandController
    %Kennedy, C. A., & Carpenter, M. H. (2016). 
    %Diagonally Implicit Runge-Kutta methods for ordinary differential equations, a review 
    
    methods
        function obj = H321(varargin)
            a = [-1/3, -1/18, 5/18];
            b = [5/6, 1/6];
            
            AQfunc = @(q, A) (A / q);
            
            obj = obj@matlode.stepsizecontroller.SoderlandController('A', a, 'B', b, 'AQfunc', AQfunc, varargin{:});
        end
    end
end

