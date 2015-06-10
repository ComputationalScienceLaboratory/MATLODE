%%
%  Authored by Tony D'Augustine, Adrian Sandu, and Hong Zhang.
%  Computational Science Laboratory, Virginia Tech.
%  ©2015 Virginia Tech Intellectual Properties, Inc.
%
function [ z ] = lss_mul_jac2( z, g, djdt )
%LSS_MUL_JAC2: DO NOT USE THIS FUNCTION. ONLY HERE TO MAKE TRANSLATING
%EASIER

    z = djdt*g;

end

