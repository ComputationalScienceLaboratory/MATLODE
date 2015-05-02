function [ err ] = relativeError( y,y_ref, p )

    if ( nargin < 3 )
        p = 2;
    end
    
    err = norm(y-y_ref,p)/norm(y_ref,p);
    
return;

