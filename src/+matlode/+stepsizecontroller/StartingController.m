classdef (Abstract) StartingController < matlode.stepsizecontroller.StepSizeController
    %Adaptive controller
    % This class is mostly used to contian the starting procedure of all
    % adaptive controllers and some properties
    
    properties (Constant)
        Adaptive = true;
    end
    
    properties (SetAccess = immutable)
        Fac
        FacMax
        FacMin
        InitalStep
    end
    
    
    methods
        function obj = StartingController(hist, varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            
            p.addParameter('Fac', 0.8);
            p.addParameter('FacMin', 0.1);
            p.addParameter('FacMax', 2);
            p.addParameter('InitalStep', matlode.startingstep.BookStarting());
            
            p.parse(varargin{:});
            
            opts = p.Results;
            varargout = p.Unmatched;
            
           obj = obj@matlode.stepsizecontroller.StepSizeController(hist, varargout); 
           
            obj.Fac = opts.Fac;
            obj.FacMin = opts.FacMin;
            obj.FacMax = opts.FacMax;
            obj.InitalStep = opts.InitalStep;
            
        end
        
        function [h0, f0, fevals] = startingStep(obj, f, tspan, y0, order, errFunc, minStep, maxStep)
            [h0, f0, fevals] = obj.InitalStep.startingStep(f, tspan, y0, order, errFunc, minStep, maxStep);
        end
    end
end

