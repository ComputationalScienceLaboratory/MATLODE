classdef ERK < matlode.rk.RungeKutta
    %Explicit Runge Kutta class
    
    properties (Constant)
        PartitionMethod = false
    end
    
    methods
        function obj = ERK(a, b, bHat, c, e, order, embeddedOrder)
            obj = obj@matlode.rk.RungeKutta(a, b, bHat, c, e, order, embeddedOrder);
        end
    end
    
    methods (Access = protected)
        function opts = matlodeSets(obj, p, varargin)
            
            %ERK specific options
            
            opts = matlodeSets@matlode.rk.RungeKutta(obj, p, varargin{:});
        end
        
        function [k, yN, err] = timeStep(obj, f, t, y, h, k, ~, prevAccept)
            persistent fsal s a b c
            if isempty(fsal)
                fsal = obj.FsalStart;
                s = obj.Stage;
                a = obj.A;
                b = obj.B;
                c = obj.C;
            end
            
            if fsal && prevAccept
                k(:, 1) = k(:, end);
            end
            
            for i = fsal:s
                z = y;
                for j = 1:i-1
                    z = z + k(:, j) .* (h .* a(i, j));
                end
                thc = t + h .* c(i);
                k(:, i) = f(thc, z);
            end
            yN = y + k * (h .* b');
            err = 1;
        end
        
        function [k, yN, err] = timeStepErr(obj,  f, t, y, h, k, ErrNorm, prevAccept)
            persistent fsal s a b c e
            if isempty(fsal)
                fsal = obj.FsalStart;
                s = obj.Stage;
                a = obj.A;
                b = obj.B;
                c = obj.C;
                e = obj.E;
            end
            
            if fsal && prevAccept
                k(:, 1) = k(:, end);
            end
            
            for i = fsal:s
                z = y;
                for j = 1:i-1
                    z = z + k(:, j) .* (h .* a(i, j));
                end
                thc = t + h .* c(i);
                k(:, i) = f(thc, z);
            end
            kh = k .* h;
            yN = y + kh * b';
            yE = kh * e';
            err = ErrNorm.errEstimate([yN, y], yE);
        end
        
        
    end
end

