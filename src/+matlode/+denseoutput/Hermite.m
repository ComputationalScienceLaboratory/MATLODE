classdef Hermite < matlode.denseoutput.DenseOutput
    %HERMITE Formula obtained from Solving ODEs I Chapter II.6
    
    methods 
        function obj = Hermite()
        end
        
        function [denseY, fEvals] = denseOut(~, f, t, tneed, y0, y1, ~, dt)
            
            theta = (tneed - t) / h;
            f0 = f.F(t, y0);
            f1 = f.F(t + h, y1);
            denseY = (1 - theta) * y0 + theta * y1 + theta*(theta - 1) * ((1-2*theta) * (y1 - y0) + (theta - 1) * (dt * f0 + theta * dt * f1));
            fEvals = 2;
        end
    end
end

