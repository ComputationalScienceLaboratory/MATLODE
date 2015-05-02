function [ Tout, Yout ] = matlOde_allocateMemory( nvar, chunk )

    Tout = zeros(chunk, 1);
    Yout = zeros(nvar, chunk);

return;

