classdef (Abstract) StepSizeController < handle
    %Template for a step size controller
    
    properties (SetAccess = immutable)
        History
        Fac
        FacMax
        FacMin
        InitalStep
    end
    
    methods (Access = protected)
        function obj = StepSizeController(hist, varargin)
            
            p = inputParser;
            p.KeepUnmatched = false;
            
            p.addParameter('Fac', 0.8, matlode.util.scalarValidationFunc);
            p.addParameter('FacMin', 0.1, matlode.util.scalarValidationFunc);
            p.addParameter('FacMax', 2, matlode.util.scalarValidationFunc);
            p.addParameter('InitalStep', matlode.startingstep.WattsStarting());
            
            %Gives better error checking
            p.parse(varargin{:});
            
            opts = p.Results;
           
            obj.Fac = opts.Fac;
            obj.FacMin = opts.FacMin;
            obj.FacMax = opts.FacMax;
            obj.InitalStep = opts.InitalStep;
            
           obj.History = hist; 
        end
    end
    
    methods (Abstract)
        [accept, hNew, tNew] = newStepSize(obj, prevAccept, t, tspan, h, err, q);
	end
    
    methods
        function [h0, f0, fevals] = startingStep(obj, f, tspan, y0, order, errFunc, minStep, maxStep)
            [h0, f0, fevals] = obj.InitalStep.startingStep(f, tspan, y0, order, errFunc, minStep, maxStep);
        end
    end
    
end

