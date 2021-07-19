classdef Hermite < matlode.denseoutput.DenseOutput
    %HERMITE Formula obtained from Solving ODEs I Chapter II.6
    
    methods 
        function obj = Hermite()
        end
        
        function [denseY, fEvals] = denseOut(~, f, t, tNeeded, yC, yN, ~, h)
            
            theta = (tNeeded - t) / h;
            f0 = f(t, yC);
            f1 = f(t + h, yN);
            denseY = (1 - theta) * yC + theta * yN + theta*(theta - 1) * ((1-2*theta) * (yN - yC) + (theta - 1) * (h * f0 + theta * h * f1));
            fEvals = 2;
        end
    end
end

