function [ Tout, Yout ] = matlOde_appendMemory( nvar, t, y, chunk )

    [ Tout, Yout ] = matlOde_allocateMemory( nvar, chunk );
    
    Tout = [ t; Tout ];
    Yout = [ y Yout ];

return;

