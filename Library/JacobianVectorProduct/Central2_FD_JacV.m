function [ JacV, Stats ] = Central2_FD_JacV( T, Y, V, OdeFunction, Stats )
    yDimension = length(Y);

    % Central finite difference (order 2)
    JacV = zeros(yDimension,1);
    canonical = eye(yDimension);
    epsilon = sqrt((1+norm(Y))*eps)/norm(V);    
    for i=1:yDimension
        JacV = JacV + ((OdeFunction(T,Y+epsilon*canonical(:,i))-OdeFunction(T,Y-epsilon*canonical(:,i)))/(2*epsilon))*V(i);
    end
    Stats.Nfun = Stats.Nfun + 2*yDimension;
    
return;