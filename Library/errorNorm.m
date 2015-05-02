function Err = errorNorm( NVAR, Y, SCAL )

% Vectorized    
    Err = sum((Y.*SCAL).^2,1);    
    Err = max( sqrt( Err/double(NVAR) ), 1.0d-10 );

    
return;

