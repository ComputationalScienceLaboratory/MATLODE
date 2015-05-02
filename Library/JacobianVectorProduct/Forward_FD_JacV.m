function [ JacV, Stats ] = Forward_FD_JacV( T, Y, V, F0, OdeFunction, Stats )
    yDimension = length(Y);

    % Forward finite difference
    JacV = zeros(yDimension,1);
    canonical = eye(yDimension);
    %epsilon = (1/(yDimension*norm(V)))*sum(sqrt(eps).*Y+sqrt(eps));    
    epsilon = sqrt((1+norm(Y))*eps)/norm(V);
    for i=1:yDimension
        JacV = JacV + ((OdeFunction(T,Y+epsilon*canonical(:,i))-F0)/epsilon)*V(i);
    end
    Stats.Nfun = Stats.Nfun + yDimension;
    
return;