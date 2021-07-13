classdef BookStarting < matlode.startingstep.StartingStep
    %This starting step procedure is derived from the Solving ODEs I book
    
    
    methods
        function obj = BookStarting(varargin)
            obj = obj@matlode.startingstep.StartingStep(varargin{:});
        end
        
        function [h0, f0, fevals] = startingStep(obj, f, tspan, y0, order, errFunc, minStep, maxStep)
            fevals = 2;
            tdir = sign(tspan(end) - tspan(1));
            
            f0 = f(tspan(1), y0);
            
            d0 = errFunc([y0, y0], 0);
            d1 = errFunc([f0, f0], 0);
            
            if any([d0, d1] < 1e-5)
                h0 = 1e-6;
            else
                h0 = 0.01 * d0 / d1;
            end

            y1 = y0 + h0 * f0 * tdir;
            f1 = f(tspan(1) + h0 * tdir, y1);
            err0 = errFunc([f1 - f0, f1 - f0], 0);
            d2 = err0 / h0;
            dm = max(d1, d2);

            if dm <= eps
                h1 = max(1e-6, h0 * 1e-3);
            else
                h1 = (0.01 / dm)^(1 / (order + 1));
            end
            
            
            h0 = min(100 * h0, h1);
            h0 = max([min([h0, maxStep * tdir]), minStep * tdir]) * tdir;
        end
    end
end

