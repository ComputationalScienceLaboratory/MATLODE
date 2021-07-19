classdef (Abstract) DenseOutput < handle
    methods (Abstract)
        [denseY, fEvals] = denseOut(obj, f, t, tNeeded, yC, yN, k, h)
    end
end

