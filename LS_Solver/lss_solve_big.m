function [ Xout ] = lss_solve_big( Transp, Zbig, ax_big, ISING )
%LSS_SOLVE_BIG Summary of this function goes here
%   Detailed explanation goes here

    Xout = ax_big\Zbig;
    Xout = Xout';


end

