function [ djdt ] = lss_jac_time_derivative( djdt, fjac, delta )
%LSS_JAC_TIME_DERIVATIVE Summary of this function goes here
%   Detailed explanation goes here

    djdt = djdt - fjac;
    a = 1/delta;
    djdt = a*djdt;

return;

