function [ JacV, Stats ] = Complex_FD_JacV( T, Y, V, OdeFunction, Stats )
    yDimension = length(Y);

    % Complex finite difference
    JacV = zeros(yDimension,1);
    canonical = eye(yDimension);
    %epsilon = (1/(yDimension*norm(V)))*sum(sqrt(eps).*Y+sqrt(eps));    
    epsilon = sqrt((1+norm(Y))*eps)/norm(V);
    for i=1:yDimension
        JacV = JacV + ((imag(OdeFunction(T,Y+1i*canonical(:,i)*epsilon))/epsilon))*V(i);
    end
    Stats.Nfun = Stats.Nfun + yDimension;
    
return;
