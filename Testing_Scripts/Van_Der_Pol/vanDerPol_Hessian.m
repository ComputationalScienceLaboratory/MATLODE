function H = vanDerPol_Hessian(~, y)

mu = 10.0;    

H = [ 0, 0; ...
      0, 0 ];
    
H(:,:,2) = [ 2*mu*y(2), -2*mu*y(1); ...
            -2*mu*y(1), 0          ];
     
return;

