%%
%  Authored by Tony D'Augustine, Adrian Sandu, and Hong Zhang.
%  Computational Science Laboratory, Virginia Tech.
%  ©2015 Virginia Tech Intellectual Properties, Inc.
%
function [ z ] = lss_mul_jac( g, fjac )
%LSS_MUL_JAC Summary of this function goes here
%   Detailed explanation goes here

    z = fjac*g;

end

