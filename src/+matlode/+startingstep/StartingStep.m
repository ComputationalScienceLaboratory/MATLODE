classdef StartingStep
    %This is an abstract setup for dealing with different kinds of starting
    %step procedures
    
    methods
        function obj = StartingStep(varargin)
            p = inputParser;
            p.KeepUnmatched = false;
            
            %error checking using the parser
            p.parse(varargin{:});
        end
    end
    
    methods (Abstract)
        [h0, f0, fevals] = startingStep(obj, f, tspan, y0, order, errFunc, minStep, maxStep);
    end
end

