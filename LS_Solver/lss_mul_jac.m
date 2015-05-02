function [ z ] = lss_mul_jac( g, fjac )
%LSS_MUL_JAC Summary of this function goes here
%   Detailed explanation goes here

    z = fjac*g;

end

