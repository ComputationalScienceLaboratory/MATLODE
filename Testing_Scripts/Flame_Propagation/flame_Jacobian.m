function [ jac ] = flame_Jacobian( t,y )
%FLAME_JACOBIAN Summary of this function goes here
%   Detailed explanation goes here

    jac = 2*y - 3*y^2;

end

