classdef SoderlandController < matlode.stepsizecontroller.StartingController
    %SOBERLANDCONTROLLER 
    %A general controller for all history based non-conditional controllers
    %can be expanded up to any amount of history gien proper coefficents
    
    %Söderlind, G. (2003). Digital filters in adaptive time-stepping.
    % ACM Transactions on Mathematical Software 
    
    properties (SetAccess = immutable)
        A
        B
        AQfunc % Geneall @(q, A) (A / q)
    end
    
    methods
        function obj = SoderlandController(varargin)
            
            p = inputParser;
            p.KeepUnmatched = true;
            
            p.addParameter('A', [1]);
            p.addParameter('B', [0]);
            p.addParameter('AQfunc', @(q, A) (A/q));
            
            p.parse(varargin{:});
            
            opts = p.Results;
            varargout = p.Unmatched;
            
            if isempty(opts.A)
                error('A cannot be empty');
            end
            
            if isempty(opts.B)
                error('B cannots be empty');
            end
            
            if ~isa(opts.AQfunc, 'function_handle')
                error('AQfunc must be a function handle')
            end
            
            obj = obj@matlode.stepsizecontroller.StartingController(max(length(opts.A), length(opts.B)), varargout);
            
            obj.A = opts.A;
            obj.B = opts.B;
            obj.AQfunc = opts.AQfunc;
        end
        
        function [accept, hNew, tNew] = newStepSize(obj, ~, t, ~, h, err, q, ~, ~)
            
            accept = err(1) <= 1; 
            if accept
                tNew = t + h(1);
                
            else
                tNew = t;
            end
            
            
            facScal = prod(err(1:obj.History).^obj.AQfunc(q,obj.A));
            hScal = prod(abs(h(1:obj.History - 1) ./ h(2:obj.History)).^obj.B);
            prediction = obj.Fac * facScal * hScal;
            
            hNew = h(1) * min(obj.FacMax, max(obj.FacMin,  prediction));
            
            
        end
    end
end

