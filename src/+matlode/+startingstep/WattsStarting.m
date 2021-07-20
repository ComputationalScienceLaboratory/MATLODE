classdef WattsStarting < matlode.startingstep.StartingStep
    %Based of Watt's starting procedure
    %Watts, H. A. (1983). Starting step size for an ODE solver. Journal of Computational and Applied Mathematics
    
    properties (SetAccess = immutable)
        InitTol
        InitEpsMax
        InitEpsMin
        Fac
    end
    
    methods
        function obj = WattsStarting(varargin)
            
            p = inputParser;
            p.KeepUnmatched = true;
            
            p.addParameter('InitTol', eps^(1/2), matlode.util.scalarValidationFunc);
            p.addParameter('InitEpsMax', eps^(1/2), matlode.util.scalarValidationFunc);
            p.addParameter('InitEpsMin', 100 * eps, matlode.util.scalarValidationFunc);
            p.addParameter('Fac', 0.8, matlode.util.scalarValidationFunc);
            
            p.parse(varargin{:});
            
            opts = p.Results;
            varargout = p.Unmatched;
           
            obj = obj@matlode.startingstep.StartingStep(varargout);
            
            obj.InitTol = opts.InitTol;
            obj.InitEpsMax = opts.InitEpsMax;
            obj.InitEpsMin = opts.InitEpsMin;
            obj.Fac = opts.Fac;
            
        end
        
        function [h0, f0, fevals] = startingStep(obj, f, tspan, y0, order, ErrNorm, minStep, maxStep)
            fevals = 2;
            tdir = sign(tspan(end) - tspan(1));
            
            f0 = f(tspan(1), y0);
            fnorm = ErrNorm.errEstimate([f0, f0], 0);
            
            hstep = abs(tspan(1) - tspan(end));
            
            h0 = max([min([obj.InitEpsMax * abs(tspan(1)), hstep]), obj.InitEpsMin * abs(tspan(1))]);
            
            if h0 < eps
                h0 = eps^obj.InitEpsMax * hstep;
            end
            
            y1 = y0 + h0 * tdir * f0;
            
            ydd = (1/(h0)) * (f(tspan(1) + h0 * tdir, y1) - f0) ;
            
            fddnorm = ErrNorm.errEstimate([ydd, ydd], 0);
            
            if fddnorm > eps
                h0 = obj.Fac * sqrt(2) * obj.InitTol^(1/(order + 1)) / sqrt(fddnorm);
            elseif fnorm > eps
                h0 = obj.Fac * obj.InitTol^(1/(order + 1)) / fnorm;
            else
                h0 = obj.Fac * obj.InitTol^(1/(order + 1)) * abs(tspan(2) - tspan(1));
            end
            
            h0 = max([min([h0, maxStep * tdir]), minStep * tdir]) * tdir;
        end
    end
end

