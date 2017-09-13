function result = phiFun(A, LA, UA, pA, dtA, gA)
    M = length(A);

    if ~exist('LA', 'var') || ~exist('UA', 'var') ...
            || ~exist('pA', 'var') || ~exist('dtA', 'var') || ~exist('gA','var')
        result = (expm(A) - eye(M))/(A);
    else
        phi1_Nr = (expm(A) - eye(M));
        result = UA\(LA\(phi1_Nr(pA, :)/(dtA * gA)));
    end
return