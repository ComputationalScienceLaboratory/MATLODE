classdef (Abstract) StartingController < matlode.stepsizecontroller.StepSizeController
    %Adaptive controller
    % This class is mostly used to contian the starting procedure of all
    % adaptive controllers and some properties
    
    methods
        function [h0, f0] = startingStep(obj, f, tspan, y0, order, errFunc, intialStep)
            if intialStep ~= 0
                f0 = f(tspan(1), y0);
                h0 = intialStep;
            else
                f0 = f(tspan(1), y0);
                d0 = errFunc([y0, y0], 0);
                d1 = errFunc([f0, f0], 0);
                if any([d0, d1] < 1e-5)
                    h0 = 1e-6;
                else
                    h0 = 0.01 * d0 / d1;
                end

                y1 = y0 + h0 * f0;
                f1 = f(tspan(1) + h0, y1);
                err0 = errFunc([f1 - f0, f1 - f0], 0);
                d2 = err0 / h0;
                dm = max(d1, d2);

                if dm <= 1e-15
                    h1 = max(1e-6, h0 * 1e-3);
                else
                    h1 = (0.01 / dm)^(1 / (order + 1));
                end

                h0 = min(100 * h0, h1);
            end
            
            %allocate for length of history
            h0 = ones(1, obj.History) * h0;
        end
    end
end

