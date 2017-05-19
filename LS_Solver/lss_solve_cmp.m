%%
%  Authored by Tony D'Augustine, Adrian Sandu, and Hong Zhang.
%  Computational Science Laboratory, Virginia Tech.
%  ©2015 Virginia Tech Intellectual Properties, Inc.
%
function [ R2_real R3_imag ] = lss_solve_cmp( NVAR, Transp, R2_real, R3_imag, e2 )
%LSS_SOLVE_CMP Summary of this function goes here

    rhs = zeros(NVAR,1);

    for m=1:NVAR
        rhs(m) = R2_real(m) + R3_imag(m)*1i;
    end
    X = e2\rhs;
    %Xout = X';
    
    R2_real = real(X);
    R3_imag = imag(X);

return;

