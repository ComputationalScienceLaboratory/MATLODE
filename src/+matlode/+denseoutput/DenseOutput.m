classdef (Abstract) DenseOutput < handle
    methods (Abstract)
        [denseY, fEvals] = denseOut(obj, f, t, tneed, y0, y1, stages, dt)
    end
end

