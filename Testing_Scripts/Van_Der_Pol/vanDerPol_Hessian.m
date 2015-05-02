function H = vanDerPol_Hessian( ~,y,mu )

    if nargin < 3 %supply default if not given
        mu = 10.0;
    end


    H = [ 0, 0; ...
          0, 0 ];
    
    H(:,:,2) = [ 2*mu*y(2), -2*mu*y(1); ...
                -2*mu*y(1), 0          ];
     
return;

