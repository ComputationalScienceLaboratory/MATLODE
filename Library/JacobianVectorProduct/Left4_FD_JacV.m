function [ JacV, Stats ] = Left4_FD_JacV( T, Y, V, OdeFunction, Stats )
    yDimension = length(Y);

    % Central finite difference (order 4)
    JacV = zeros(yDimension,1);
    canonical = eye(yDimension);
    epsilon = sqrt((1+norm(Y))*eps)/norm(V);    
    for i=1:yDimension
        JacV = JacV + ((-OdeFunction(T,Y+3*epsilon*canonical(:,i))-6*OdeFunction(T,Y+2*epsilon*canonical(:,i)) ...
            +18*OdeFunction(T,Y+epsilon*canonical(:,i))-10*OdeFunction(T,Y) ...
            -3*OdeFunction(T,Y-epsilon*canonical(:,i)))/(12*epsilon))*V(i);
    end
    Stats.Nfun = Stats.Nfun + 4*yDimension;
    
return;