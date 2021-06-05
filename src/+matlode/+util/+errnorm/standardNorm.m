function err = standardNorm(y, yHat, options)
        %procedure described in solving ODEs I Harrier
        %cannot assume that user has signal processing toolbox

        sc = (options.AbsTol + max(abs(y), abs(yHat)) .* options.RelTol);
        
        err = sqrt(1/length(y)) * norm((y - yHat) ./ sc);
end

