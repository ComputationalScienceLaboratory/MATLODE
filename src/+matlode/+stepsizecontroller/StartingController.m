classdef (Abstract) StartingController < matlode.stepsizecontroller.StepSizeController
    % This class is used to contain common adaptive stepsizecontroler
    % properties and to allow the user to change the starting step
    % procedure
    
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
            
            p.addParameter('Fac', 0.8, matlode.util.scalarValidationFunc);
            p.addParameter('FacMin', 0.1, matlode.util.scalarValidationFunc);
            p.addParameter('FacMax', 2, matlode.util.scalarValidationFunc);
            p.addParameter('InitalStep', matlode.startingstep.WattsStarting());
            
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

