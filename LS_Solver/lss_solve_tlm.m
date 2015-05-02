function [ ISING Xout ] = lss_solve_tlm( trans, rhs, e_tlm, ISING )
%LSS_SOLVE_TLM Summary of this function goes here
%   Detailed explanation goes here

    X = rhs\e_tlm;
    Xout = X';

end

