classdef SoberlandController < matlode.stepsizecontroller.StartingController
    %SOBERLANDCONTROLLER 
    %A general controller for all history based non-conditional controllers
    %can be expanded up to any amount of history gien proper coefficents
    
    properties (Constant)
        Adaptive = true;
    end
    
    properties (SetAccess = immutable)
        Fac
        FacMax
        FacMin
        A
        B
        AQfunc % Geneall @(q, A) (A / (1+q))
        History
    end
    
    methods
        function obj = SoberlandController(varargin)
            
             p = inputParser;
            
            p.addParameter('Fac', 0.95);
            p.addParameter('FacMin', 0.1);
            p.addParameter('FacMax', 2);
            p.addParameter('B', [0]);
            p.addParameter('A', [1]);
            p.addRequired('AQfunc');
            
            p.parse(varargin{:});
            
            opts = p.Results;
            
            obj.Fac = opts.Fac;
            obj.FacMin = opts.FacMin;
            obj.FacMax = opts.FacMax;
            
            if isempty(opts.A)
                error('A cannot be empty')
            end
            
            obj.A = opts.A;
            obj.B = opts.B;
            
            %This is need because input parser struggles with
            %function_handles
            if ~isa(opts.AQfunc, 'function_handle')
                error('AQfunc must be a function handle')
            end
            
            obj.AQfunc = opts.AQfunc;
            obj.History = length(obj.A);
        end
        
        function [accept, hNew, tNew] = newStepSize(obj, ~, t, ~, h, err, q, ~, ~)
            
            accept = mean(err) <= 1; 
            if accept
                tNew = t + h(1);
                
            else
                tNew = t;
            end
            
            
            facScal = prod(err(1:obj.History).^obj.AQfunc(q,obj.A));
            hScal = prod((h(1:obj.History - 1) ./ h(2:obj.History)).^obj.B);
            prediction = obj.Fac * facScal * hScal;
            
            hNew = h(1) * min(obj.FacMax, max(obj.FacMin,  prediction));
            
            
        end
    end
end

