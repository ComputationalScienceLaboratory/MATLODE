classdef Gustafson < matlode.stepsizecontroller.StartingController
    %Gustafsson, K. (1991). Control theoretic techniques for stepsize 
    % selection in explicit Runge-Kutta methods. ACM Transactions on Mathematical Software
    
    
    properties (SetAccess = immutable)
        Ki
        Kp
        A
    end
    
    methods
        function obj = Gustafson(varargin)
            
            p = inputParser;
            p.KeepUnmatched = true;
            
            p.addParameter('Ki', 0.3, matlode.util.scalarValidationFunc);
            p.addParameter('Kp', -0.4, matlode.util.scalarValidationFunc);
            p.addParameter('A', 1, matlode.util.scalarValidationFunc);
            
            p.parse(varargin{:});
            
            opts = p.Results;
            varargout = p.Unmatched;
            
            obj = obj@matlode.stepsizecontroller.StartingController(2, varargout);
            
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

