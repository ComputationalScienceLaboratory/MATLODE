classdef SetStep < matlode.startingstep.StartingStep
    %SETSTEP Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        StartingStep
    end
    
    methods
        function obj = SetStep(ssp, varargin)
            obj = obj@matlode.startingstep.StartingStep(varargin{:});
            obj.StartingStep = ssp;
        end
        
        function [h0, f0, fevals] = startingStep(obj, ~, tspan, ~, ~, ~, minStep, maxStep)
            fevals = 0;
            f0 = [];
            if sign(obj.StartingStep) ~= sign(tspan(end) - tspan(1))
                error('The starting step cannot go in the opposite direction, %f', obj.StartingStep)
            end
            if obj.StartingStep > maxStep
               error('The starting step is greater then the maximum step %f', maxStep) 
            elseif obj.StartingStep < minStep
                error('The starting step cannot be smaller than the minimum step %f', minStep)
            end
            
            h0 = max([min([obj.StartingStep, maxStep]), minStep]);
        end
    end
end

