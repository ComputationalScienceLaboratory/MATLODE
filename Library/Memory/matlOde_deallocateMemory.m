function [ Tout, Yout ] = matlOde_deallocateMemory( Tout, Yout, Nacc )

    Tout(Nacc+1:end,:) = [];
    Yout(:,Nacc+1:end) = [];

return;

