%%
%  Authored by Tony D'Augustine, Adrian Sandu, and Hong Zhang.
%  Computational Science Laboratory, Virginia Tech.
%  ©2015 Virginia Tech Intellectual Properties, Inc.
%
function [ z ] = lss_mul_jac1( z, g, fjac1 )
%LSS_MUL_JAC1: DO NOT USE THIS FUNCTION. ONLY HERE TO MAKE TRANSLATING
%EASIER
%   Detailed explanation goes here

    z = fjac1*g;

end

