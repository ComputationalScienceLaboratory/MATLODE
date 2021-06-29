classdef Gustafson < matlode.stepsizecontroller.StartingController
    %TODO Reference Gustaffson paper(1991)
    %Slight derivation of the Gustafson paper on the coefficents
    
    properties (Constant)
        Adaptive = true;
    end
    
    properties (SetAccess = immutable)
        Fac
        FacMin
        FacMax
        Ki
        Kp
        A
        History
    end
    
    methods
        function obj = Gustafson(varargin)
            
            obj.History = 2;
            
            p = inputParser;
            
            p.addParameter('Fac', 0.95);
            p.addParameter('FacMin', 0.1);
            p.addParameter('FacMax', 2);
            p.addParameter('Ki', 0.3);
            p.addParameter('Kp', -0.4);
            p.addParameter('A', 1);
            
            p.parse(varargin{:});
            
            opts = p.Results;
            
            obj.Fac = opts.Fac;
            obj.FacMin = opts.FacMin;
            obj.FacMax = opts.FacMax;
            obj.Ki = opts.Ki;
            obj.Kp = opts.Kp;
            obj.A = opts.A;
            
        end
        
        function [accept, hNew, tNew] = newStepSize(obj, prevAccept, t, ~, h, err, q, ~, ~)
            accept = err(1) <= 1;
            
            scalh = 1;
            
            if accept
                tNew = t + h(1);
                if ~prevAccept
                    %H_0220 Controller
                   scalh = (h(1) / h(2))^obj.A; 
                end
                %PI Controller
                scalerr = err(1)^(obj.Ki/q) * err(2)^(obj.Kp/q);
            else
                tNew = t;
                %Standard Controller(I controller)
                scalerr = err(1)^(-1/(q + 1));
            end
            
            hNew = h(1) * min(obj.FacMax, max(obj.FacMin, obj.Fac * scalerr * scalh));
            
        end
    end
end

